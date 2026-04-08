#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(matrixStats)
})

######################################################################
# Helpers
######################################################################

fread_auto <- function(path) fread(path, na.strings=c("NA","NaN","","-9"))
trim_split <- function(x) trimws(unlist(strsplit(x, ",")))

# brace expansion: PC{1:10} -> PC1..PC10
expand_braces <- function(tokens, avail) {
  out <- character(0)
  for (tk in tokens) {
    m <- regexec("^([A-Za-z0-9_.]+)\\{([0-9]+):([0-9]+)\\}$", tk)
    g <- regmatches(tk, m)[[1]]
    if (length(g)==4) {
      base <- g[2]; a <- as.integer(g[3]); b <- as.integer(g[4])
      seqn <- if (a<=b) a:b else b:a
      out <- c(out, paste0(base, seqn))
    } else out <- c(out, tk)
  }
  miss <- setdiff(out, avail)
  if (length(miss))
    stop("Covariate(s) not found: ", paste(miss, collapse=", "))
  unique(out)
}

######################################################################
# Options
######################################################################

option_list <- list(
  make_option("--bfile", type="character",
              help="PRE-FILTERED PLINK prefix (.bed/.bim/.fam)"),
  make_option("--covar", type="character", default=NULL,
              help="Covariate file (FID IID + columns)"),
  make_option("--covar_cont_cols", type="character", default=NULL),
  make_option("--covar_bin_cols",  type="character", default=NULL),
  make_option("--pheno", type="character", help="Phenotype file"),
  make_option("--pheno_cols", type="character",
              help="Phenotypes to define Step-1 sample"),
  make_option("--trait", type="character", default="qt", help="qt|bt"),

  # for bt: union-of-fails across all foci
  make_option("--bt_weight_cols", type="character", default=NULL,
              help="Binary phenotypes for weighting; default=--pheno_cols"),

  make_option("--bsize", type="integer", default=1000),
  make_option("--numtol", type="double", default=1e-6,
              help="Variance gate threshold"),
  make_option("--r2_covar_max", type="double", default=0.999999),
  make_option("--probe_snp", type="character", default=NULL),
  make_option("--threads", type="integer",
              default=parallel::detectCores(logical=FALSE)),
  make_option("--out_prefix", type="character",
              default="regenie_step1_exact_R12")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("\n========================================================\n")
cat("  REGENIE Step-1 SNP Filter (R12)\n")
cat("  * No internal keep/extract — use PLINK pre-filtered bfile\n")
cat("  * Exact variance gate: VAR(post-projection) < numtol\n")
cat("  * Binary: union of fails across multiple BT phenotypes\n")
cat("========================================================\n\n")

######################################################################
# Load BIM/FAM
######################################################################

bim <- fread(paste0(opt$bfile,".bim"),
  col.names=c("CHR","SNP","CM","BP","A1","A2"))
fam <- fread(paste0(opt$bfile,".fam"),
  col.names=c("FID","IID","PID","MID","SEX","PHENO"))
fam[,KEY := paste(FID,IID,sep="_")]

######################################################################
# Phenotypes & sample
######################################################################

stopifnot(!is.null(opt$pheno_cols))
ph_list <- trim_split(opt$pheno_cols)

ph <- fread_auto(opt$pheno)
stopifnot(all(c("FID","IID") %in% names(ph)))
stopifnot(all(ph_list %in% names(ph)))

ph[,KEY := paste(FID,IID,sep="_")]
ph[,OK := apply(.SD, 1, function(x) all(!is.na(x))), .SDcols=ph_list]
ph_ok <- ph[OK==TRUE]

sample_dt <- merge(
  fam[,.(FID,IID,KEY)],
  ph_ok[,.(KEY)],
  by="KEY", all=FALSE
)

if (!nrow(sample_dt))
  stop("No samples remain after phenotype intersection.")

######################################################################
# Covariates
######################################################################

has_cov <- !is.null(opt$covar) &&
           (!is.null(opt$covar_cont_cols) || !is.null(opt$covar_bin_cols))

cont_cols <- bin_cols <- character(0)

if (has_cov) {
  cv <- fread_auto(opt$covar)
  stopifnot(all(c("FID","IID") %in% names(cv)))
  cv[,KEY := paste(FID,IID,sep="_")]

  cname <- setdiff(names(cv), c("FID","IID","KEY"))

  cont_cols <- if (!is.null(opt$covar_cont_cols)) trim_split(opt$covar_cont_cols) else character(0)
  bin_cols  <- if (!is.null(opt$covar_bin_cols))  trim_split(opt$covar_bin_cols)  else character(0)

  if (length(cont_cols)) cont_cols <- expand_braces(cont_cols, cname)
  if (length(bin_cols))  bin_cols  <- expand_braces(bin_cols, cname)

  if (length(intersect(cont_cols, bin_cols)))
    stop("Covariate appears in both continuous and categorical lists.")

  used <- unique(c(cont_cols, bin_cols))

  sample_dt <- merge(sample_dt,
                     cv[,c("KEY",used),with=FALSE],
                     by="KEY", all=FALSE)

  # drop NA covariates
  na_rows <- sample_dt[,apply(.SD,1,anyNA),.SDcols=used]
  sample_dt <- sample_dt[!na_rows]
}

if (!nrow(sample_dt))
  stop("No samples remain after phenotype/covariate merge.")

ind_row <- match(sample_dt$KEY, fam$KEY)
if (anyNA(ind_row))
  stop("Internal error: NA in ind_row; FAM mismatch.")

N <- length(ind_row)
cat("Step-1 sample size =", N, "\n")

######################################################################
# Design matrix C
######################################################################

build_C <- function() {
  n <- nrow(sample_dt)
  if (!has_cov) return(cbind(Intercept=rep(1,n)))

  Xc <- NULL
  if (length(cont_cols)) {
    tmp <- as.data.frame(sample_dt[,..cont_cols])
    if (!all(sapply(tmp,is.numeric)))
      stop("Continuous covariates must be numeric.")
    Xc <- as.matrix(tmp)
    colnames(Xc) <- cont_cols
  }

  Xb <- NULL
  if (length(bin_cols)) {
    tmpb <- sample_dt[,..bin_cols]
    tmpb[] <- lapply(tmpb,
        function(v){ if(is.logical(v)) v<-as.integer(v); as.factor(v) })
    MM <- model.matrix(~., data=as.data.frame(tmpb))
    MM <- MM[,-1,drop=FALSE]
    Xb <- as.matrix(MM)
  }

  X0 <- cbind(Xc, Xb)
  if (is.null(X0)) return(cbind(Intercept=rep(1,n)))

  C <- cbind(Intercept=rep(1,n), X0)
  qrC <- qr(C)
  C[,qrC$pivot[seq_len(qrC$rank)],drop=FALSE]
}

C <- build_C()
stopifnot(nrow(C) == N)

######################################################################
# Binary or quantitative weighting plans
######################################################################

if (opt$trait=="bt") {

  # which BT phenotypes to weight on?
  bt_list <- if (!is.null(opt$bt_weight_cols))
               trim_split(opt$bt_weight_cols) else ph_list

  if (!all(bt_list %in% names(ph)))
    stop("Some BT weight phenotypes not in phenotype file.")

  if (length(bt_list)<1)
    stop("Trait=bt but no binary phenotypes available.")

  # prepare:
  QW_list <- vector("list", length(bt_list)); names(QW_list) <- bt_list
  SW_list <- vector("list", length(bt_list)); names(SW_list) <- bt_list
  fail_by_bt <- vector("list", length(bt_list)); names(fail_by_bt) <- bt_list

  # fit null logistic per BT phenotype
  for (k in seq_along(bt_list)) {
    y <- ph_ok[match(sample_dt$KEY,ph_ok$KEY), get(bt_list[k])]
    if (!all(y %in% c(0,1)))
      stop("BT phenotypes must be 0/1; offending=", bt_list[k])

    gl <- suppressWarnings(glm(y ~ C[,-1,drop=FALSE], family=binomial))
    mu <- fitted(gl)
    W  <- pmax(mu*(1-mu), .Machine$double.eps)
    sw <- sqrt(W)

    Cw <- C * sw
    qrCw <- qr(Cw)
    Qw <- qr.Q(qr(Cw[,qrCw$pivot[seq_len(qrCw$rank)],drop=FALSE]))

    QW_list[[k]] <- Qw
    SW_list[[k]] <- sw
  }

} else {
  Q_unw <- qr.Q(qr(C))
}

######################################################################
# Load FBM from pre-filtered BED
######################################################################

rds <- paste0(opt$out_prefix,".rds")
if (!file.exists(rds)) {
  cat("Converting BED -> FBM …\n")
  snp_readBed2(paste0(opt$bfile,".bed"), backingfile=opt$out_prefix)
}
obj <- snp_attach(rds)
G <- obj$genotypes

if (ncol(G)!=nrow(bim))
  stop("FBM columns != BIM rows (check pre-filtered bfile).")

if (max(ind_row)>nrow(G))
  stop("ind_row exceeds FBM rows (bfile mismatch).")

######################################################################
# Step-1 blockwise variance processing
######################################################################

candidate <- seq_len(nrow(bim))
bsize <- opt$bsize
blocks <- split(candidate, ceiling(seq_along(candidate)/bsize))

fail_union <- integer(0)

probe_done <- FALSE
blk_i <- 0L

for (idx in blocks) {
  blk_i <- blk_i + 1

  if (!length(idx)) next

  # load genotypes
  X <- as.matrix(G[ind_row, idx, drop=FALSE])

  # allowed codes check
  if (any(!(X %in% c(0,1,2,3) | is.na(X)))) {
    bad <- X[!(X %in% c(0,1,2,3) | is.na(X))][1]
    stop("Illegal genotype code: ", bad)
  }

  # mean-impute code 3
  miss <- (X==3 & !is.na(X))
  if (any(miss)) {
    X2 <- X; X2[miss] <- NA_real_
    m  <- colMeans(X2, na.rm=TRUE)
    for (j in seq_len(ncol(X)))
      if (any(miss[,j])) X[miss[,j],j] <- m[j]
  }

  # impute any NA
  if (anyNA(X)) {
    m2 <- colMeans(X, na.rm=TRUE)
    for (j in seq_len(ncol(X))) {
      nas <- is.na(X[,j])
      if (any(nas)) X[nas,j] <- m2[j]
    }
  }

  # --- binary traits: union of fails ---
  if (opt$trait=="bt") {

    for (k in seq_along(bt_list)) {
      Qw <- QW_list[[k]]
      sw <- SW_list[[k]]

      Xw  <- X * sw
      Xcw <- sweep(Xw, 2, colMeans2(Xw), "-")

      if (nrow(Qw) != nrow(Xcw))
        stop("Internal: nrow(Qw)!=nrow(Xcw) for ", bt_list[k])

      Proj <- Qw %*% (crossprod(Qw, Xcw))
      Rw   <- Xcw - Proj
      Rw   <- sweep(Rw, 2, colMeans2(Rw), "-")

      var_res <- colVars(Rw)
      failed  <- which(is.na(var_res) | (var_res < opt$numtol))

      if (length(failed)) {
        fail_union <- c(fail_union, idx[failed])
        fail_by_bt[[bt_list[k]]] <-
          c(fail_by_bt[[bt_list[k]]], idx[failed])
      }

      # per-BT probe
      if (!is.null(opt$probe_snp)) {
        hit <- which(bim$SNP[idx] == opt$probe_snp)
        if (length(hit)) {
          j <- hit[1]
          cat(sprintf("PROBE %s [%s]: col=%d var_res=%.10g decision=%s\n",
                      opt$probe_snp, bt_list[k], idx[j], var_res[j],
                      if (is.na(var_res[j])||var_res[j]<opt$numtol) "FAIL" else "KEEP"))
        }
      }
    }

  } else {

    # --- quantitative ---
    Xcw <- sweep(X, 2, colMeans2(X), "-")

    if (nrow(Q_unw) != nrow(Xcw))
      stop("Internal: nrow(Q_unw)!=nrow(Xcw).")

    Proj <- Q_unw %*% (crossprod(Q_unw, Xcw))
    Rw   <- Xcw - Proj
    Rw   <- sweep(Rw, 2, colMeans2(Rw), "-")

    var_res <- colVars(Rw)
    failed  <- which(is.na(var_res) | (var_res < opt$numtol))
    if (length(failed)) fail_union <- c(fail_union, idx[failed])

    # QT probe
    if (!probe_done && !is.null(opt$probe_snp)) {
      hit <- which(bim$SNP[idx] == opt$probe_snp)
      if (length(hit)) {
        j <- hit[1]
        cat(sprintf("PROBE %s: col=%d var_res=%.10g decision=%s\n",
                    opt$probe_snp, idx[j], var_res[j],
                    if (is.na(var_res[j])||var_res[j]<opt$numtol) "FAIL" else "KEEP"))
        probe_done <- TRUE
      }
    }

  }

  if (blk_i %% 50 == 0)
    cat(sprintf(" ... processed block %d (up to SNP index %d)\n",
                blk_i, max(idx)))

}

######################################################################
# Final keep list (union of fails removed)
######################################################################

fail_union <- sort(unique(fail_union))
keep_idx  <- setdiff(seq_len(nrow(bim)), fail_union)

######################################################################
# Output
######################################################################

dir.create(dirname(opt$out_prefix), showWarnings=FALSE, recursive=TRUE)
f_keep <- paste0(opt$out_prefix,".snps_keep.txt")
f_low  <- paste0(opt$out_prefix,".snps_lowvar_union.txt")
f_sum  <- paste0(opt$out_prefix,".summary.tsv")

fwrite(data.table(SNP=bim$SNP[keep_idx]), f_keep, col.names=FALSE)
fwrite(data.table(SNP=bim$SNP[fail_union]), f_low, col.names=FALSE)

# per-BT diagnostic
if (opt$trait=="bt") {
  for (nm in names(fail_by_bt)) {
    if (!length(fail_by_bt[[nm]])) next
    out_bt <- paste0(opt$out_prefix,".snps_lowvar_",nm,".txt")
    fwrite(data.table(SNP=bim$SNP[sort(unique(fail_by_bt[[nm]]))]),
           out_bt, col.names=FALSE)
  }
}

fwrite(data.table(
  metric=c("N_total","N_fail_union","N_kept","numtol","trait",
           if(opt$trait=="bt") "N_BT_traits" else NULL),
  value=c(nrow(bim), length(fail_union), length(keep_idx),
          opt$numtol, opt$trait,
          if(opt$trait=="bt") length(bt_list) else NULL)
), f_sum, sep="\t")

cat("\n========================================================\n")
cat(" DONE.\n")
cat(" Keep (safe for all binary traits): ", f_keep,"\n")
cat(" Low-var (union):                   ", f_low,"\n")
if (opt$trait=="bt")
  cat(" Per-BT low-var lists written.\n")
cat(" Summary:                           ", f_sum,"\n")
cat("========================================================\n\n")
