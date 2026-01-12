version 1.1

import "regenie.wdl"

workflow group_pheno {

  input { 
    Array[Array[File]] genos
    File? covar
    String? covarColList
    String? catCovarList
    String? phenoColList
    File phenotype_table
    File phewas_manifest
    Float missing_thresh
  }

  String prefix = basename(phenotype_table, ".tsv.gz")

  scatter(geno in genos) {
    call regenie.filter_genos {
      input:
        bed = geno[0],
        bim = geno[1],
        fam = geno[2],
        pheno = phenotype_table
    }
  }

  call regenie.merge_genos {
    input:
      beds = filter_genos.out_bed,
      bims = filter_genos.out_bim,
      fams = filter_genos.out_fam,
      prefix = prefix
  }

  call regenie.filter_snps as filter_snps {
    input:
      bed = merge_genos.out_bed,
      bim = merge_genos.out_bim,
      fam = merge_genos.out_fam,
      prefix = prefix
  }

  call split_phenotypes {
    input:
      phenotype_table = phenotype_table,
      phewas_manifest = phewas_manifest,
      covar = covar,
      phenoColList = phenoColList
  }

  call calc_dist as calc_dist_bt {
    input:
      pheno = split_phenotypes.bin,
      qc_id = filter_snps.qc_id
  }
  
  call cluster_traits as cluster_traits_bt {
    input:
      pheno = split_phenotypes.bin,
      dist = calc_dist_bt.out,
      missing_thresh = missing_thresh,
      qc_id = filter_snps.qc_id
  }
  
  call group_ids as group_ids_bt {
    input:
      pheno = split_phenotypes.bin,
      qc_id = filter_snps.qc_id,
      pheno_groups = cluster_traits_bt.pheno_groups,
      covar = covar
  }
  
  scatter(g in zip(cluster_traits_bt.pheno_groups, group_ids_bt.group_ids)) {
    String group_id_bt = basename(g.right, ".ids")

    call regenie.filter_snps as filter_snps_bt {
      input:
        bed = merge_genos.out_bed,
        bim = merge_genos.out_bim,
        fam = merge_genos.out_fam,
        samples_keep = g.right,
        prefix = group_id_bt
    }

    call regenie.step1 as step1_bt {
      input:
        bed = merge_genos.out_bed,
        bim = merge_genos.out_bim,
        fam = merge_genos.out_fam,
        pheno = split_phenotypes.bin,
        phenoColList = g.left,
        covar = covar,
        qc_id = g.right,
        qc_snplist = filter_snps_bt.qc_snplist,
        covarColList = covarColList,
        catCovarList = catCovarList,
        bt = true,
        prefix = group_id_bt
    }
  }

  call calc_dist as calc_dist_qt {
    input:
      pheno = split_phenotypes.quant,
      qc_id = filter_snps.qc_id
  }
  
  call cluster_traits as cluster_traits_qt {
    input:
      pheno = split_phenotypes.quant,
      dist = calc_dist_qt.out,
      missing_thresh = missing_thresh,
      qc_id = filter_snps.qc_id
  }

  call group_ids as group_ids_qt {
    input:
      pheno = split_phenotypes.quant,
      qc_id = filter_snps.qc_id,
      pheno_groups = cluster_traits_qt.pheno_groups,
      covar = covar
  }
  
  scatter(g in zip(cluster_traits_qt.pheno_groups, group_ids_qt.group_ids)) {
    String group_id_qt = basename(g.right, ".ids")

    call regenie.filter_snps as filter_snps_qt {
      input:
        bed = merge_genos.out_bed,
        bim = merge_genos.out_bim,
        fam = merge_genos.out_fam,
        samples_keep = g.right,
        prefix = group_id_qt
    }

    call regenie.step1 as step1_qt {
      input:
        bed = merge_genos.out_bed,
        bim = merge_genos.out_bim,
        fam = merge_genos.out_fam,
        pheno = split_phenotypes.quant,
        phenoColList = g.left,
        covar = covar,
        qc_id = g.right,
        qc_snplist = filter_snps_qt.qc_snplist,
        covarColList = covarColList,
        catCovarList = catCovarList,
        bt = false,
        prefix = group_id_qt
    }
  }

  call merge_pred_list as merge_pred_list_bt {
    input:
      pred_lists = step1_bt.pred_list
  }

  call merge_pred_list as merge_pred_list_qt {
    input:
      pred_lists = step1_qt.pred_list
  }

  output {
    Array[File] pred_list = [merge_pred_list_bt.out, merge_pred_list_qt.out]
    Array[File] loco_qt = flatten(step1_qt.loco)
    Array[File] loco_bt = flatten(step1_bt.loco)
    File pheno_bt = split_phenotypes.bin
    File pheno_qt = split_phenotypes.quant
  }
}

task split_phenotypes {

  input {
    File phenotype_table
    File phewas_manifest
    File? covar
    String? phenoColList
  }

  String out_bin = basename(phenotype_table, ".gz") + "_bin.txt"
  String out_quant = sub(out_bin, "bin", "quant")

  command <<<
    Rscript - <<-SPLIT_PHENO
      library(tidyverse)
      manifest <- read_csv("~{phewas_manifest}",
                          col_types = cols(field_code="d",
                                      QC_flag_ID="d",
                                          date_code="d",
                                            age_code="d",
                                            included_in_analysis="d",
                                            quant_combination="d",
                                            lower_limit="d",
                                            upper_limit="d",
                                            .default="c"))
      id_bin <- manifest %>%
        filter(analysis == "binary") %>%
        pull(PheWAS_ID)

      id_quant <- manifest %>%
        filter(analysis %in% c("quant", "count")) %>%
        pull(PheWAS_ID)

      pheno <- read_tsv("~{phenotype_table}", col_types = cols(.default = "d"))

      phenoColList <- "~{select_first([phenoColList, ''])}"

      if (phenoColList != "") {
        pheno_cols <- str_split_1(phenoColList, ",")
        pheno <- pheno %>%
          select(eid, any_of(pheno_cols))
      }

      split_pheno <- function(ids, out, covar) {
        pheno_split <- pheno %>%
          mutate(FID=eid) %>%
          select(FID, IID=eid, any_of(ids))

        if (!missing(covar)) {
          cov <- read_tsv("~{covar}", col_types = cols(.default = "d")) %>%
            drop_na %>%
            select(FID, IID)

          pheno_keep <- pheno_split %>% 
            inner_join(cov, by=c("FID", "IID")) %>%
            summarise(across(c(-FID, -IID), ~sum(. == 1, na.rm = T))) %>%
            pivot_longer(everything(), names_to = "pheno", values_to = "n_cases") %>%
            filter(n_cases >= 10) %>%
            pull(pheno)
    
          pheno_split <- pheno_split %>%
            select(FID, IID, any_of(pheno_keep))
        }

        pheno_keep <- pheno_split %>%
            summarise(across(c(-FID, -IID), ~sum(!is.na(.)))) %>%
            pivot_longer(everything(), names_to = "pheno", values_to = "n") %>%
            filter(n >= 100) %>%
            pull(pheno)
        
        pheno_split %>%
          select(FID, IID, any_of(pheno_keep)) %>%
          write_tsv(out)
      }

      split_pheno(id_bin, "~{out_bin}" ~{', covar=' + '"' + covar + '"'})
      split_pheno(id_quant, "~{out_quant}") 
    SPLIT_PHENO
  >>>

  output {
    File bin = out_bin
    File quant = out_quant
  }

  runtime {
    memory: "200 GB"
    container: "rocker/tidyverse"
  }
}

task calc_dist {

  input {
    File pheno
    File qc_id
  }

  String dist_out = basename(pheno, ".txt") + "_dist.RDS"

  command <<<
    Rscript - <<-CALC_DIST
      install.packages("parallelDist")
      library(parallelDist)
      library(tidyverse)
      
      pheno <- read_tsv("~{pheno}")
      qc_id <- read_tsv("~{qc_id}", col_names=c("FID", "IID"))

      pheno.m <- pheno %>%
        inner_join(qc_id, by=c("FID", "IID")) %>%
        select(-FID, -IID) %>%
        mutate(across(everything(), ~is.na(.) %>% ifelse(0, 1))) %>%
        as.matrix %>%
        t

      pheno_dist <- parallelDist(pheno.m, "binary", threads=12)
      saveRDS(pheno_dist, "~{dist_out}")
    CALC_DIST
  >>>

  runtime {
    memory: "200 GB"
    container: "rocker/tidyverse"
  }

  output {
    File out = dist_out
  }
}

task cluster_traits {

  input {
    File pheno
    File qc_id
    File dist
    Float missing_thresh
  }

  String groups_out = basename(pheno, ".txt") + "_group"

  command <<<
    Rscript - <<-CLUSTER_TRAITS
      install.packages("usedist")
      library(tidyverse)
      
      qc_id <- read_tsv("~{qc_id}", col_names=c("FID", "IID"))

      pheno <- read_tsv("~{pheno}") %>%
        inner_join(qc_id, by=c("FID", "IID"))

      # Keep phenotype columns in a matrix; store names
      pheno_cols <- pheno %>% select(-FID, -IID)
      phenos <- colnames(pheno_cols)

      # NA matrix (TRUE = missing, FALSE = observed)
      M <- is.na(as.matrix(pheno_cols))   # n_samples x n_phenos

      pheno_dist <- readRDS("~{dist}")
      pheno_clust <- hclust(pheno_dist, method="ward.D2")

      low <- 1
      high <- ncol(M)
      best_k <- NA

      assign_all_k <- cutree(pheno_clust, k = seq_len(high)) # precompute all clusterings
  
      missing_props_for_group <- function(idx, M) {
        rows_keep <- rowSums(!M[, idx, drop = FALSE]) > 0
        if (!any(rows_keep)) {
          rep(NA_real_, length(idx))
        } else {
          colMeans(M[rows_keep, idx, drop = FALSE])
        }
      }
  
      while (low <= high) {
        k <- floor((low + high) / 2)
        groups_k <- assign_all_k[, k]
        split_idx <- split(seq_len(ncol(M)), groups_k)
   
        cluster_missing <- lapply(split_idx, function(idx) missing_props_for_group(idx, M))
        all_ok <- all(unlist(cluster_missing) < ~{missing_thresh}, na.rm = TRUE)
   
        if (all_ok) {
          best_k <- k
          high <- k - 1  # try smaller k
        } else {
          low <- k + 1   # need more clusters
        }
      }
  
      # Return best_k and its cluster assignments
      if (!is.na(best_k)) {
        message("Optimum number of clusters = ", best_k)
        groups_k <- assign_all_k[, best_k]
        final_clusters <- tibble(pheno = phenos, group = groups_k)
      } else {
        message("No clustering satisfies the threshold.")
        final_clusters <- NULL
      }

      final_clusters %>%
        group_by(group) %>%
        summarise(pheno=paste(pheno, collapse = ",")) %>%
        ungroup %>%
        mutate(group=as.integer(group)) %>%
        arrange(group) %>%
        pull(pheno) %>%
        cat(file="~{groups_out}.txt", sep="\n")
    CLUSTER_TRAITS
  >>>

  runtime {
    memory: "64 GB"
    container: "rocker/tidyverse"
  }

  output {
    Array[String] pheno_groups = read_lines(groups_out + ".txt")
  }
}

task group_ids {

  input {
    File pheno
    File qc_id
    Array[String] pheno_groups
    File? covar
  }

  String groups_out = basename(pheno, ".txt") + "_group"
  Int n_groups = length(pheno_groups)
  File samples_keep = select_first([covar, pheno])

  command <<<
    Rscript - <<-GROUP_IDS
      library(tidyverse)
      
      qc_id <- read_tsv("~{qc_id}", col_names=c("FID", "IID"))

      pheno <- read_tsv("~{pheno}") %>%
        inner_join(qc_id, by=c("FID", "IID"))

      samples_keep <- read_tsv("~{samples_keep}") %>%
        select(FID, IID)

      pheno <- pheno %>%
        inner_join(samples_keep, by=c("FID", "IID"))

      pheno_groups <- "~{sep=';' pheno_groups}"

      final_clusters <- tibble(pheno=str_split_1(pheno_groups, ";")) %>%
        mutate(group=1:n()) %>%
        separate_longer_delim(pheno, ",")

      final_samples <- final_clusters %>%
        group_by(group) %>%
        group_map(~pheno %>%
          select(c(FID, IID, .x\$pheno)) %>%
            filter(if_any(.x\$pheno, ~!is.na(.))))

      walk(1:length(final_samples),
          ~write_tsv(final_samples[[.]] %>%
            select(FID, IID),
              paste0("~{groups_out}_", . - 1, ".ids"), col_names = FALSE))
    GROUP_IDS
  >>>

  runtime {
    memory: "64 GB"
    container: "rocker/tidyverse"
  }

  output {
    Array[File] group_ids = suffix(".ids", prefix("~{groups_out}_", range(n_groups)))
  }
}

task merge_pred_list {

  input {
    Array[File] pred_lists
  }

  String pred_list_merged = sub(basename(select_first(pred_lists)), "_group_[0-9]+", "")

  command <<<
    cat ~{sep=' ' pred_lists} > ~{pred_list_merged}
  >>>

  runtime {
    memory: "4 GB"
  }

  output {
    File out = pred_list_merged
  }
}
