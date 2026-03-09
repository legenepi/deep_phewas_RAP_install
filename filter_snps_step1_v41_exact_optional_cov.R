#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(matrixStats)
})

# ============================ Options ============================
option_list <- list(
  make_option("--bfile",              type="character", help="PLINK prefix (no .bed/.bim/.fam)"),
  make_option("--covar",              type="character", default=NULL,
              help="(Optional) Covariate file with FID IID + covariates"),
  make_option("--covar_cont_cols",    type="character", default=NULL,
              help="(Optional) Comma-separated CONTINUOUS covariates (e.g., AGE,PC1,PC2 or PC{1:10})"),
  make_option("--covar_bin_cols",     type="character", default=NULL,
              help="(Optional) Comma-separated BINARY/CATEGORICAL covariates (e.g., SEX,BATCH)"),
  make_option("--pheno",              type="character", help="Phenotype file with FID IID + traits"),
  make_option("--pheno_cols",         type="character", help="Comma-separated phenotype names (Step-1 sample = intersection)"),
  make_option("--pheno_col",          type="character", default=NULL,
              help="(Binary only) phenotype column (0/1) used to fit null logistic weights"),
  make_option("--trait",              type="character", default="qt", help="Trait type for Step-1: qt|bt"),
  make_option("--keep",               type="character", default=NULL,
              help="(Optional) 2-column file FID IID listing samples to KEEP"),
  make_option("--extract",            type="character", default=NULL,
              help="(Optional) 1-column file listing SNP IDs to INCLUDE"),
  make_option("--bsize",              type="integer",  default=1000, help="Block size of SNPs (Step-1 default: 1000)"),
  make_option("--geno_max_miss",      type="double",   default=0.10, help="Max SNP missingness before checks"),
  make_option("--numtol",             type="double",   default=1e-6, help="Step-1 'Uh-oh' variance threshold (v4.1: 1e-6)"),
  make_option("--r2_covar_max",       type="double",   default=0.999999,
              help="(Diagnostic) flag if R^2(G~C) exceeds this; ignored if no covariates beyond intercept"),
  make_option("--threads",            type="integer",  default=parallel::detectCores(logical=FALSE), help="Threads"),
  make_option("--out_prefix",         type="character", default="regenie_step1_v41_exact_optcov", help="Output prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("\n=============================================================\n")
cat("  REGENIE v4.1 Step-1 Low-Variance Filter (Block-wise; cont/bin covariates optional)\n")
cat("  + --keep (samples) and --extract (variants) support\n")
cat("=============================================================\n")
cat("bsize       :", opt$bsize, "\n")
cat("numtol      :", format(opt$numtol, digits=8), "\n")
cat("trait       :", opt$trait, "\n")
cat("miss thresh :", opt$geno_max_miss, "\n")
cat("threads     :", opt$threads, "\n\n")

# ============================ Helpers ============================
fread_auto <- function(path) fread(path, na.strings=c("NA","NaN","","-9"))
trim_split <- function(x) trimws(unlist(strsplit(x, ",")))

# Minimal brace expansion like PC{1:10} -> PC1,...,PC10 (single brace pair)
expand_braces <- function(tokens, available_cols) {
  out <- character(0)
  for (tk in tokens) {
    m <- regexec("^([A-Za-z0-9_.]+)\\{([0-9]+):([0-9]+)\\}$", tk)
    g <- regmatches(tk, m)[[1]]
    if (length(g) == 4) {
      base <- g[2]; a <- as.integer(g[3]); b <- as.integer(g[4])
      seqn <- if (a <= b) a:b else a:b
      cand <- paste0(base, seqn)
      out  <- c(out, cand)
    } else {
      out <- c(out, tk)
    }
  }
  miss <- setdiff(out, available_cols)
  if (length(miss)) stop("Covariate column(s) not found after brace expansion: ", paste(miss, collapse=", "))
  unique(out)
}

# Read --keep (FID IID), with/without header
read_keep_keys <- function(path) {
  dt <- fread(path, header = FALSE, col.names = paste0("V", seq_len(max(count.fields(path, sep="\t")))))
  # try to infer columns
  nm <- names(dt)
  if (all(c("FID","IID") %in% nm)) {
    k <- paste(dt$FID, dt$IID, sep="_")
  } else {
    if (ncol(dt) < 2) stop("--keep requires 2 columns: FID IID")
    k <- paste(dt[[1]], dt[[2]], sep="_")
  }
  unique(k)
}

# Read --extract (SNP IDs), with/without header
read_extract_ids <- function(path) {
  dt <- fread(path, header = FALSE)
  v <- unique(dt[[1]])
  as.character(v)
}

# ============================ Load PLINK maps ============================
bim <- fread(paste0(opt$bfile, ".bim"),
             col.names=c("CHR","SNP","CM","BP","A1","A2"))
fam <- fread(paste0(opt$bfile, ".fam"),
             col.names=c("FID","IID","PID","MID","SEX","PHENO"))
fam[, KEY := paste(FID, IID, sep="_")]

# ---------- (Optional) --keep: pre-filter fam ----------
if (!is.null(opt$keep)) {
  keep_keys <- read_keep_keys(opt$keep)
  fam <- fam[KEY %in% keep_keys]
  if (nrow(fam) == 0) stop("After --keep, no individuals remain in FAM.")
}

# ---------- (Optional) --extract: predefine variant whitelist ----------
ind_col_whitelist <- seq_len(nrow(bim))
if (!is.null(opt$extract)) {
  extract_snps <- read_extract_ids(opt$extract)
  ind_col_whitelist <- which(bim$SNP %in% extract_snps)
  if (!length(ind_col_whitelist)) stop("None of the --extract SNPs were found in BIM.")
  cat("Initial SNPs in BIM:", nrow(bim), " | After --extract:", length(ind_col_whitelist), "\n")
}

# ============================ Phenotypes ============================
stopifnot(!is.null(opt$pheno), !is.null(opt$pheno_cols))
ph_cols <- trim_split(opt$pheno_cols)

pheno <- fread_auto(opt$pheno)
stopifnot(all(c("FID","IID", ph_cols) %in% names(pheno)))
pheno[, KEY := paste(FID, IID, sep="_")]
pheno[, COMPLETE := apply(.SD, 1, function(x) all(!is.na(x))), .SDcols=ph_cols]
pheno <- pheno[COMPLETE == TRUE]

if (opt$trait == "bt") {
  stopifnot(!is.null(opt$pheno_col))
  stopifnot(opt$pheno_col %in% names(pheno))
}

# ============================ Covariates (optional; cont vs bin) ============================
has_covars <- !is.null(opt$covar) && (!is.null(opt$covar_cont_cols) || !is.null(opt$covar_bin_cols))
cont_cols <- character(0); bin_cols <- character(0)

if (has_covars) {
  covar <- fread_auto(opt$covar)
  stopifnot(all(c("FID","IID") %in% names(covar)))
  covar[, KEY := paste(FID, IID, sep="_")]

  cont_cols <- if (!is.null(opt$covar_cont_cols)) trim_split(opt$covar_cont_cols) else character(0)
  bin_cols  <- if (!is.null(opt$covar_bin_cols))  trim_split(opt$covar_bin_cols)  else character(0)

  # brace expansion separately
  all_cols <- setdiff(names(covar), c("FID","IID","KEY"))
  if (length(cont_cols)) cont_cols <- expand_braces(cont_cols, all_cols)
  if (length(bin_cols))  bin_cols  <- expand_braces(bin_cols,  all_cols)

  if (length(intersect(cont_cols, bin_cols)) > 0)
    stop("A covariate cannot be both continuous and binary: ",
         paste(intersect(cont_cols, bin_cols), collapse=", "))

  used_cov_cols <- c(cont_cols, bin_cols)

  # Build Step-1 sample: fam (possibly --keep-filtered) ∩ covars ∩ phenotype-complete
  sample_dt <- merge(
    merge(fam[, .(FID,IID,KEY)],
          covar[, c("KEY", used_cov_cols), with=FALSE],
          by="KEY", all=FALSE),
    pheno[, .(KEY)],
    by="KEY", all=FALSE
  )
  # drop missing covariates
  na_row <- sample_dt[, apply(.SD, 1, function(x) any(is.na(x))), .SDcols=used_cov_cols]
  sample_dt <- sample_dt[!na_row]

} else {
  # no covariates → sample = fam (possibly --keep-filtered) ∩ phenotype-complete
  sample_dt <- merge(fam[, .(FID,IID,KEY)],
                     pheno[, .(KEY)],
                     by="KEY", all=FALSE)
}

stopifnot(nrow(sample_dt) > 0)
ind_row <- match(sample_dt$KEY, fam$KEY)
N <- length(ind_row)
cat("Unified Step-1 sample size:", N, "\n")
cat("Covariates present? ", has_covars, "\n\n")

# ============================ Design matrix C (intercept + [cont + one-hot(bin)]) ============================
build_C_matrix <- function() {
  if (!has_covars) return(cbind(Intercept = rep(1, N)))

  # continuous
  Xc <- NULL
  if (length(cont_cols)) {
    tmp <- as.data.frame(sample_dt[, ..cont_cols])
    bad <- names(tmp)[!sapply(tmp, is.numeric)]
    if (length(bad)) stop("Continuous covariates not numeric: ", paste(bad, collapse=", "))
    Xc <- as.matrix(tmp); colnames(Xc) <- cont_cols
  }
  # binary/categorical → treatment-coded dummies
  Xb <- NULL
  if (length(bin_cols)) {
    tmpb <- sample_dt[, ..bin_cols]
    tmpb[] <- lapply(tmpb, function(v) { if (is.logical(v)) v <- as.integer(v); as.factor(v) })
    MM <- model.matrix(~ . , data = as.data.frame(tmpb))
    MM <- MM[, -1, drop=FALSE]  # drop intercept
    Xb <- as.matrix(MM)
  }
  X0 <- cbind(Xc, Xb)
  if (is.null(X0)) return(cbind(Intercept = rep(1, N)))
  C <- cbind(Intercept = rep(1, N), X0)
  qrC <- qr(C)
  C[, qrC$pivot[seq_len(qrC$rank)], drop=FALSE]
}
C <- build_C_matrix()

# ============================ Binary-trait weighting (null logistic) ============================
if (opt$trait == "bt") {
  y <- pheno[match(sample_dt$KEY, pheno$KEY), get(opt$pheno_col)]
  stopifnot(all(y %in% c(0,1)))
  gl <- suppressWarnings(glm(y ~ C[,-1,drop=FALSE], family=binomial(link="logit")))
  mu <- as.numeric(fitted(gl))
  W  <- pmax(mu * (1 - mu), .Machine$double.eps)
  sw <- sqrt(W)
  Cw <- C * sw
  qrCw <- qr(Cw)
  Qw   <- qr.Q(qr(Cw[, qrCw$pivot[seq_len(qrCw$rank)], drop=FALSE]))
} else {
  Qw <- qr.Q(qr(C))
}

# ============================ Read BED as FBM ============================
rds_path <- paste0(opt$out_prefix, ".rds")
if (!file.exists(rds_path)) {
  cat("Converting BED → FBM (one-time)…\n")
  snp_readBed2(paste0(opt$bfile, ".bed"), backingfile=opt$out_prefix)
}
obj <- snp_attach(rds_path)
G <- obj$genotypes
stopifnot(ncol(G) == nrow(bim))

# ============================ Missingness screen (only on --extract whitelist) ============================
candidate_cols <- ind_col_whitelist
cat("Candidate SNPs before missingness screen:", length(candidate_cols), "\n")

cat("Computing SNP missingness…\n")
miss_count <- rep(0L, ncol(G))
for (idx in split(candidate_cols, ceiling(seq_along(candidate_cols)/2000))) {
  X <- as.matrix(G[ind_row, idx, drop=FALSE])
  miss <- (X == 3)
  miss_count[idx] <- colSums(miss)
}
miss_rate <- miss_count / N
keep_by_miss <- rep(FALSE, ncol(G))
keep_by_miss[candidate_cols] <- (miss_rate[candidate_cols] <= opt$geno_max_miss)
ind_col_all <- which(keep_by_miss)
cat("SNPs removed for missingness >", opt$geno_max_miss, ": ", length(candidate_cols) - length(ind_col_all), "\n\n")

# ============================ Exact block-wise residualization & scaling ============================
cat("Running block-wise residualization & scaling…\n")
bsize <- opt$bsize
blocks <- split(ind_col_all, ceiling(seq_along(ind_col_all) / bsize))

fail_lowvar <- logical(ncol(G))
fail_r2     <- logical(ncol(G))   # diagnostic (ignored if no covariates beyond intercept)

blk_i <- 0L
for (idx in blocks) {
  blk_i <- blk_i + 1L
  X <- as.matrix(G[ind_row, idx, drop=FALSE])

  # mean-impute (3 -> mean)
  miss <- (X == 3)
  if (any(miss)) {
    X2 <- X; X2[miss] <- NA_real_
    m  <- colMeans2(X2, na.rm=TRUE)
    for (j in seq_len(ncol(X))) if (any(miss[,j])) X[miss[,j], j] <- m[j]
  }

  # (bt only) weight the block
  Xw <- if (opt$trait == "bt") X * sw else X

  # center each SNP (block-wise)
  Xcw <- sweep(Xw, 2, colMeans2(Xw), "-")

  # project out covariates (or intercept-only; Qw handles both)
  Proj <- Qw %*% (crossprod(Qw, Xcw))
  Rw   <- Xcw - Proj

  # re-center residuals (block-wise) and compute SD
  Rw <- sweep(Rw, 2, colMeans2(Rw), "-")
  sd_resid <- sqrt(pmax(colVars(Rw), 0))  # sample SD

  # Step-1 "Uh-oh" gate
  fail_lowvar[idx] <- (sd_resid < opt$numtol) | is.na(sd_resid)

  # Diagnostic R^2(G~C) only if covariates exist beyond intercept
  has_covariate_cols <- (ncol(C) > 1)
  if (has_covariate_cols) {
    numer <- colSums(Proj^2)
    denom <- colSums(Xcw^2) + .Machine$double.eps
    r2    <- pmin(pmax(numer / denom, 0), 1)
    fail_r2[idx] <- r2 > opt$r2_covar_max
  }

  if ((blk_i %% 50) == 0)
    cat(sprintf(" ... processed %d blocks (%d SNPs)\n", blk_i, max(idx)))
}

keep_mask <- keep_by_miss & !fail_lowvar
snps_fail_lowvar <- bim$SNP[which(fail_lowvar)]
snps_keep         <- bim$SNP[which(keep_mask)]
snps_fail_r2      <- if (ncol(C) > 1) bim$SNP[which(fail_r2 & keep_by_miss)] else character(0)

# ============================ Outputs ============================
dir.create(dirname(opt$out_prefix), showWarnings=FALSE, recursive=TRUE)
f_excl <- paste0(opt$out_prefix, ".snps_excluded_lowvar.txt")
f_keep <- paste0(opt$out_prefix, ".snps_keep.txt")
f_diag <- paste0(opt$out_prefix, ".snps_flagged_r2.txt")
f_sum  <- paste0(opt$out_prefix, ".summary.tsv")

fwrite(data.table(SNP = snps_fail_lowvar), f_excl, col.names=FALSE)
fwrite(data.table(SNP = snps_keep), f_keep, col.names=FALSE)
if (ncol(C) > 1) fwrite(data.table(SNP = snps_fail_r2), f_diag, col.names=FALSE)

fwrite(data.table(
  metric=c("N_total_in_BIM","N_candidates_after_extract","N_removed_miss","N_fail_lowvar",
           if (ncol(C) > 1) "N_flagged_r2_diag" else NULL,
           "N_kept","numtol","bsize","N_individuals","N_covariates","trait","has_covariates"),
  value=c(nrow(bim), length(ind_col_whitelist), length(ind_col_whitelist) - length(ind_col_all),
          length(snps_fail_lowvar),
          if (ncol(C) > 1) length(snps_fail_r2) else NULL,
          length(snps_keep), opt$numtol, opt$bsize, N, ncol(C), opt$trait, ncol(C) > 1)
), f_sum, sep="\t")

cat("\n=============================================================\n")
cat("DONE\n")
cat("Excluded (Uh‑oh): ", f_excl, "\n")
cat("Keep list      : ", f_keep, "\n")
if (ncol(C) > 1) cat("Diag (high R^2): ", f_diag, "\n")
cat("Summary        : ", f_sum, "\n")
cat("Use with REGENIE step 1:  --extract ", f_keep, "\n", sep="")
cat("=============================================================\n\n")
