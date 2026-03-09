batch_by_comissingness <- function(
    pheno,
    meta_csv,
    prefix,
    max_size = 30,
    logfile = "batch_by_comissigness.log"
){
  
  # ====================== Packages ======================
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(Matrix)
    library(cluster)
    library(readr)
    library(purrr)
    library(ggplot2)
  })
  
  # ====================== Logging =======================
  write("", logfile)  # reset
  log_msg <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(...))
    cat(msg, "\n")
    write(msg, file=logfile, append=TRUE)
  }
  
  log_msg("==== Starting batch_by_comissingness v3.3.1 (patched) ====")
  
  # ====================== 1) Load phenotype TSV ======================
  log_msg("Loading phenotype")
  Y_raw <- pheno
  sample_id <- Y_raw[[1]]
  Y <- as.data.frame(Y_raw[,-1])
  rownames(Y) <- sample_id
  pheno_ids <- colnames(Y)
  log_msg("Loaded ", nrow(Y), " samples × ", ncol(Y), " phenotypes.")
  
  # ====================== 2) Load metadata ===========================
  log_msg("Loading metadata: ", meta_csv)
  meta_raw <- readr::read_csv(meta_csv, col_types = cols(field_code="d",
                                                         QC_flag_ID="d",
                                                         date_code="d",
                                                         age_code="d",
                                                         included_in_analysis="d",
                                                         quant_combination="d",
                                                         lower_limit="d",
                                                         upper_limit="d",
                                                         .default="c"))
  
  meta <- meta_raw %>%
    rename(
      phenotype_id = PheWAS_ID,
      domain       = pheno_group,
      analysis_type = analysis
    )
  
  if ("included_in_analysis" %in% names(meta)) {
    n0 <- nrow(meta)
    meta <- meta %>% filter(included_in_analysis == 1)
    log_msg("Filtered metadata: ", n0, " → ", nrow(meta))
  }
  
  # ====================== 3) _age mapping ============================
  log_msg("Mapping _age phenotypes to base IDs...")
  base_names <- gsub("_age$", "", pheno_ids)
  original_colnames <- pheno_ids
  
  meta <- meta %>% filter(phenotype_id %in% base_names)
  
  meta_ordered <- map_df(base_names, function(b) {
    idx <- which(meta$phenotype_id == b)
    if (length(idx) == 0) stop("Missing metadata for phenotype: ", b)
    meta[idx[1], ]
  })
  
  meta <- meta_ordered
  meta$original_colname <- original_colnames
  meta$base_name <- base_names
  
  # >>> PATCH: add unique, stable column position for each phenotype column
  meta$col_pos <- seq_len(nrow(meta))
  # <<< END PATCH
  
  # ====================== 4) Type ================================
  meta$type <- ifelse(meta$analysis_type %in% c("quant","count"),
                      "quantitative", "binary")
  log_msg("Type distribution: ",
          paste(capture.output(print(table(meta$type))), collapse=" "))
  
  # ====================== 5) Global missingness + case fraction ====
  log_msg("Computing global missingness & case_fraction ...")
  meta$missingness <- colMeans(is.na(Y))
  meta$case_fraction <- NA_real_
  for (j in seq_len(ncol(Y))) {
    if (meta$type[j] == "binary") {
      x <- Y[[j]]; x_non_na <- x[!is.na(x)]
      if (length(unique(x_non_na)) <= 2 && all(x_non_na %in% c(0,1))) {
        meta$case_fraction[j] <- mean(x_non_na == 1)
      } else if ("case_code" %in% names(meta)) {
        cc <- meta$case_code[j]
        if (!is.na(cc)) meta$case_fraction[j] <- mean(x == cc, na.rm=TRUE)
      }
    }
  }
  
  png(paste0(prefix, "missingness_hist_unrestricted.png"), width=1200, height=800)
  print(
    ggplot(meta, aes(missingness)) +
      geom_histogram(bins=50, fill="#3182ce") +
      theme_minimal(base_size=14) +
      labs(title="Global Missingness (Unrestricted)",
           x="Missing fraction", y="Phenotype count")
  )
  dev.off()
  
  # ====================== 6) Strata ================================
  assign_quant_stratum <- function(m){
    if (m < 0.05) "Q1" else if (m < 0.10) "Q2" else if (m < 0.15) "Q3" else "Q4"
  }
  assign_binary_stratum <- function(cf){
    if (is.na(cf)) NA else if (cf < 0.01) "B4" else if (cf < 0.05) "B3" else if (cf < 0.15) "B2" else "B1"
  }
  meta$quant_stratum  <- sapply(meta$missingness, assign_quant_stratum)
  meta$binary_stratum <- sapply(meta$case_fraction, assign_binary_stratum)
  
  # ====================== 7) Missingness binning ===================
  log_msg("Assigning missingness bins ...")
  miss_bins <- c(0.00, 0.15, 0.40, 1.00)
  miss_bin_labels <- c("LOW","MID","HIGH")
  bin_max_missing <- c(LOW=0.15, MID=0.30, HIGH=0.60)
  
  assign_miss_bin <- function(m) {
    cut(m, breaks=miss_bins, labels=miss_bin_labels, include_lowest=TRUE, right=TRUE)
  }
  meta$miss_bin <- assign_miss_bin(meta$missingness)
  log_msg("Missingness bin distribution: ",
          paste(capture.output(print(table(meta$miss_bin))), collapse=" "))
  
  # ====================== 8) Jaccard co-missingness ================
  log_msg("Computing Jaccard co-missingness matrix ...")
  Ms <- Matrix(is.na(as.matrix(Y))*1L, sparse=TRUE)
  colnames(Ms) <- base_names
  C  <- crossprod(Ms)
  Cd <- as.matrix(C)
  mv <- diag(Cd)
  Umat <- outer(mv, mv, "+") - Cd
  Jmat <- ifelse(Umat > 0, Cd/Umat, 0)
  diag(Jmat) <- 1
  D_all_mat <- 1 - Jmat
  
  png(paste0(prefix, "_jaccard_heatmap_unrestricted.png"), width=1600, height=1600)
  stats::heatmap(Jmat, symm=TRUE, main="Jaccard Co-Missingness (Unrestricted)")
  dev.off()
  
  # ====================== 9) Build subgroups (domain×type×strata×bin) =====
  log_msg("Building subgroups ...")
  subgroups <- list()
  for (bin in miss_bin_labels) {
    # Quantitative
    for (dom in unique(meta$domain)) {
      for (q in c("Q1","Q2","Q3","Q4")) {
        g <- meta %>% filter(miss_bin==bin, domain==dom, type=="quantitative", quant_stratum==q)
        if (nrow(g) > 0) { g$miss_bin <- bin; subgroups[[length(subgroups)+1]] <- g }
      }
    }
    # Binary
    for (dom in unique(meta$domain)) {
      for (b in c("B1","B2","B3","B4")) {
        g <- meta %>% filter(miss_bin==bin, domain==dom, type=="binary", binary_stratum==b)
        if (nrow(g) > 0) { g$miss_bin <- bin; subgroups[[length(subgroups)+1]] <- g }
      }
    }
  }
  log_msg("Total subgroups: ", length(subgroups))
  
  # ====================== 10) Clustering engine (safe PAM + recursion) =====
  
  pam_safe <- function(D_sub, k) {
    n <- attr(D_sub, "Size")
    if (is.null(n) || n < 2) return(NULL)
    if (k < 1 || k >= n) return(NULL)
    out <- try(cluster::pam(D_sub, k=k, diss=TRUE), silent=TRUE)
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  }
  
  force_split_large_cluster <- function(cl_df, max_size) {
    n <- nrow(cl_df); if (n <= max_size) return(list(cl_df))
    idx <- split(seq_len(n), ceiling(seq_len(n)/max_size))
    lapply(idx, function(i) cl_df[i,,drop=FALSE])
  }
  
  cluster_sample_indices <- function(Y, cl_df) {
    pcs <- cl_df$original_colname
    keep <- rowSums(!is.na(Y[, pcs, drop=FALSE])) > 0
    which(keep)
  }
  
  restricted_stats <- function(Y, cl_df) {
    pcs <- cl_df$original_colname
    Sidx <- cluster_sample_indices(Y, cl_df)
    if (length(Sidx) == 0) {
      cl_df$missingness_restricted <- NA_real_
      stats_df <- data.frame(n_pheno=nrow(cl_df), n_samples=0,
                             mean_miss_restr=NA_real_,
                             median_miss_restr=NA_real_,
                             max_miss_restr=NA_real_)
      return(list(cl=cl_df, stats=stats_df, samples=character(0)))
    }
    miss_restr <- colMeans(is.na(Y[Sidx, pcs, drop=FALSE]))
    cl_df$missingness_restricted <- miss_restr
    stats_df <- data.frame(n_pheno=nrow(cl_df), n_samples=length(Sidx),
                           mean_miss_restr=mean(miss_restr),
                           median_miss_restr=median(miss_restr),
                           max_miss_restr=max(miss_restr))
    return(list(cl=cl_df, stats=stats_df, samples=rownames(Y)[Sidx]))
  }
  
  pam_split_indices <- function(idxs, D_all_mat, k) {
    if (length(idxs) < 2) return(list(idxs))
    D_sub <- as.dist(D_all_mat[idxs, idxs, drop=FALSE])
    pam_fit <- pam_safe(D_sub, k)
    if (is.null(pam_fit)) return(list(idxs))
    split(idxs, pam_fit$clustering)
  }
  
  refine_cluster_recursive <- function(idxs, meta, Y, D_all_mat,
                                       max_size, target_missing,
                                       depth=0, max_depth=20,
                                       prev_mean_missing=NULL) {
    if (depth > max_depth) return(list(meta[idxs,,drop=FALSE]))
    cl_df <- meta[idxs,,drop=FALSE]
    parts <- force_split_large_cluster(cl_df, max_size)
    out <- list()
    for (part in parts) {
      rs <- restricted_stats(Y, part)
      part <- rs$cl; stats <- rs$stats; mean_r <- stats$mean_miss_restr
      if (is.na(mean_r) || mean_r <= target_missing || nrow(part)==1) {
        out <- c(out, list(part)); next
      }
      if (!is.null(prev_mean_missing) && abs(mean_r - prev_mean_missing) < 1e-6) {
        out <- c(out, list(part)); next
      }
      # >>> PATCH: use the unique column positions instead of matching base names
      part_idxs <- part$col_pos
      # <<< END PATCH
      splits <- pam_split_indices(part_idxs, D_all_mat, k=2)
      if (length(splits)==1 || any(sapply(splits, length) == nrow(part))) {
        out <- c(out, list(part)); next
      }
      for (sp in splits) {
        out <- c(out, refine_cluster_recursive(sp, meta, Y, D_all_mat,
                                               max_size, target_missing,
                                               depth+1, max_depth, mean_r))
      }
    }
    out
  }
  
  pam_constrained_by_size <- function(g_df, base_names_master, D_all_mat, max_size) {
    if (nrow(g_df) <= max_size) return(list(g_df))
    # >>> PATCH: use the unique column positions rather than matching by base name
    idx <- g_df$col_pos
    # <<< END PATCH
    low <- 1; high <- min(nrow(g_df), max(1, ceiling(nrow(g_df)/10))+10); best <- NULL
    while (low <= high) {
      k <- floor((low + high)/2)
      if (k >= nrow(g_df)) { low <- k + 1; next }
      splits <- pam_split_indices(idx, D_all_mat, k)
      if (all(sapply(splits, length) <= max_size)) { best <- splits; high <- k - 1 }
      else { low <- k + 1 }
    }
    if (is.null(best)) return(list(g_df))
    # Map split indices (col positions) back to rows in g_df
    lapply(best, function(ix) g_df[match(ix, g_df$col_pos), , drop=FALSE])
  }
  
  # ====================== 11) Run clustering across subgroups =========
  log_msg("Clustering subgroups with bin-specific targets ...")
  all_clusters <- list(); cid <- 1
  for (sg_i in seq_along(subgroups)) {
    g <- subgroups[[sg_i]]
    bin_lab <- unique(g$miss_bin); if (length(bin_lab)!=1) bin_lab <- as.character(g$miss_bin[1])
    target_missing <- bin_max_missing[[bin_lab]]
    log_msg(sprintf("Subgroup %d | bin=%s | size=%d | domain=%s | type=%s | target=%.2f",
                    sg_i, bin_lab, nrow(g), g$domain[1], g$type[1], target_missing))
    init_parts <- pam_constrained_by_size(g, base_names, D_all_mat, max_size)
    for (part in init_parts) {
      # >>> PATCH: seed refinement with the unique column positions
      p_idxs <- part$col_pos
      # <<< END PATCH
      refined <- refine_cluster_recursive(p_idxs, meta, Y, D_all_mat, max_size, target_missing)
      for (rp in refined) {
        rp$cluster_id <- cid; rp$miss_bin <- bin_lab
        all_clusters[[cid]] <- rp
        log_msg("  → Final cluster ", cid, " size=", nrow(rp), " bin=", bin_lab)
        cid <- cid + 1
      }
    }
  }
  log_msg("Clusters before outlier splitting: ", length(all_clusters))
  
  # ====================== 12) Outlier peeling ========================
  log_msg("Outlier splitting ...")
  detect_restricted_outliers <- function(cl_df, trait_max_missing=0.50, iqr_mult=1.5) {
    miss <- cl_df$missingness_restricted
    if (all(!is.finite(miss))) return(integer(0))
    abs_out <- which(miss > trait_max_missing & is.finite(miss))
    m_ok <- miss[is.finite(miss)]; Q <- quantile(m_ok, c(0.25,0.75)); IQRv <- Q[2]-Q[1]
    cutoff <- Q[2] + iqr_mult*IQRv
    iqr_out <- which(miss > cutoff & is.finite(miss))
    sort(unique(c(abs_out, iqr_out)))
  }
  split_outlier_phenotypes <- function(Y, cl_df, trait_max_missing=0.50, iqr_mult=1.5,
                                       make_singletons=TRUE, max_size=30,
                                       min_samples_per_trait=NULL, log_fn=log_msg) {
    rs0 <- restricted_stats(Y, cl_df); cl_df <- rs0$cl
    out_idx <- detect_restricted_outliers(cl_df, trait_max_missing, iqr_mult)
    if (length(out_idx)==0) return(list(clusters=list(cl_df), dropped=NULL))
    core_idx <- setdiff(seq_len(nrow(cl_df)), out_idx)
    new_clusters <- list(); dropped <- NULL
    if (length(core_idx) > 0) {
      core <- cl_df[core_idx,,drop=FALSE]
      for (part in force_split_large_cluster(core, max_size)) new_clusters <- c(new_clusters, list(part))
    }
    if (make_singletons) {
      for (j in out_idx) {
        one <- cl_df[j,,drop=FALSE]; rs1 <- restricted_stats(Y, one); nS <- rs1$stats$n_samples
        if (!is.null(min_samples_per_trait) && is.finite(nS) && nS < min_samples_per_trait) {
          dropped <- rbind(dropped, data.frame(
            phenotype = one$original_colname, base_name = one$base_name,
            n_samples = nS, missingness_restricted = rs1$stats$mean_miss_restr,
            stringsAsFactors = FALSE))
          log_fn(sprintf("Dropping outlier phenotype %s (samples %d < min %d)",
                         one$original_colname, nS, min_samples_per_trait))
        } else {
          new_clusters <- c(new_clusters, list(rs1$cl))
        }
      }
    } else {
      outliers <- cl_df[out_idx,,drop=FALSE]
      for (part in force_split_large_cluster(outliers, max_size)) new_clusters <- c(new_clusters, list(part))
    }
    list(clusters=new_clusters, dropped=dropped)
  }
  postprocess_split_outliers <- function(Y, clusters, trait_max_missing=0.50, iqr_mult=1.5,
                                         make_singletons=TRUE, max_size=30,
                                         min_samples_per_trait=NULL, log_fn=log_msg) {
    out <- list(); dropped_all <- NULL
    for (i in seq_along(clusters)) {
      res <- split_outlier_phenotypes(Y, clusters[[i]], trait_max_missing, iqr_mult,
                                      make_singletons, max_size, min_samples_per_trait, log_fn)
      out <- c(out, res$clusters)
      if (!is.null(res$dropped)) dropped_all <- rbind(dropped_all, res$dropped)
    }
    list(clusters=out, dropped=dropped_all)
  }
  post <- postprocess_split_outliers(Y, all_clusters,
                                     trait_max_missing=0.50, iqr_mult=1.5,
                                     make_singletons=TRUE, max_size=max_size,
                                     min_samples_per_trait=NULL, log_fn=log_msg)
  all_clusters <- post$clusters
  log_msg("Clusters after outlier splitting: ", length(all_clusters))
  if (!is.null(post$dropped)) {
    write.csv(post$dropped, paste0(prefix, "_dropped_outlier_phenotypes.csv"), row.names=FALSE)
    log_msg("Dropped ", nrow(post$dropped), " high-missingness phenotypes with low N.")
  }
  
  # ====================== 13) Merged missingness summary =============
  log_msg("Computing merged missingness summary ...")
  merged_stats <- list(); all_restr_miss <- c()
  for (i in seq_along(all_clusters)) {
    cl <- all_clusters[[i]]
    un_mean <- mean(cl$missingness); un_med <- median(cl$missingness); un_max <- max(cl$missingness)
    rs <- restricted_stats(Y, cl); cl <- rs$cl; all_clusters[[i]] <- cl
    re_mean <- rs$stats$mean_miss_restr; re_med <- rs$stats$median_miss_restr; re_max <- rs$stats$max_miss_restr
    nS <- rs$stats$n_samples
    all_restr_miss <- c(all_restr_miss, cl$missingness_restricted)
    merged_stats[[i]] <- data.frame(
      cluster_id = i, miss_bin = cl$miss_bin[1],
      n_pheno = nrow(cl), n_samples = nS,
      mean_missing_unrestricted = un_mean,
      median_missing_unrestricted = un_med,
      max_missing_unrestricted = un_max,
      mean_missing_restricted = re_mean,
      median_missing_restricted = re_med,
      max_missing_restricted = re_max,
      stringsAsFactors = FALSE
    )
  }
  merged_stats_df <- dplyr::bind_rows(merged_stats)
  write.csv(merged_stats_df, paste0(prefix, "_cluster_missingness_summary.csv"), row.names=FALSE)
  
  # Diagnostics: restricted missingness histogram & size distribution
  png(paste0(prefix, "_missingness_hist_restricted.png"), width=1200, height=800)
  print(
    ggplot(data.frame(miss=all_restr_miss[is.finite(all_restr_miss)]), aes(miss)) +
      geom_histogram(bins=50, fill="#2f855a") +
      theme_minimal(base_size=14) +
      labs(title="Restricted Missingness", x="Missing fraction", y="Count")
  ); dev.off()
  
  png(paste0(prefix, "_cluster_size_distribution.png"), width=1200, height=800)
  print(
    ggplot(data.frame(size=sapply(all_clusters, nrow)), aes(size)) +
      geom_histogram(binwidth=1, fill="#b83280") +
      theme_minimal(base_size=14) +
      labs(title="Final Cluster Sizes", x="# phenotypes", y="Cluster count")
  ); dev.off()
  
  # ====================== 14) Write batch & keep files ===============
  log_msg("Writing batch phenotype lists and --keep files")
  for (i in seq_along(all_clusters)) {
    cl <- all_clusters[[i]]
    writeLines(cl$original_colname, file.path(paste0("batch_", i, ".txt")))
    Sidx <- cluster_sample_indices(Y, cl); keep_ids <- rownames(Y)[Sidx]
    keep_df <- data.frame(FID=keep_ids, IID=keep_ids, stringsAsFactors=FALSE)
    write.table(keep_df, file=file.path(paste0("keep_batch_", i, ".txt")),
                col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  
  log_msg("==== Done (v3.3.1 patched) ====")
  return(all_clusters)
}
