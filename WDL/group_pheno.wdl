version 1.1

import "regenie.wdl"
import "csv_overlap.wdl" as csv

workflow group_pheno {

  input { 
    Array[Array[File]] genos
    File? covar
    String? covarColList
    String? catCovarList
    String? phenoColList
    File phenotype_table
    File phewas_manifest
    File R_batch_function
    File R_filter_lowvar_function
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

  call cluster_traits {
    input:
      pheno = phenotype_table,
      manifest = phewas_manifest,
      qc_id = filter_snps.qc_id,
      R_batch_function = R_batch_function,
      prefix = prefix
  }
  
  call split_phenotypes {
    input:
      phenotype_table = phenotype_table,
      phewas_manifest = phewas_manifest,
      covar = covar,
      phenoColList = phenoColList
  }

  scatter(g in zip(cluster_traits.batches_bt, cluster_traits.samples_bt)) {
    String group_id_bt = basename(g.left, ".txt")
    String pheno_list_bt = sep(",", read_lines(g.left))

    if (defined(phenoColList)) {
      call csv.csv_overlap as overlap_bt {
        input:
          query_csv = pheno_list_bt,
          haystack_csv = select_first([phenoColList, ""])
      }
      
      String pheno_overlap_bt = overlap_bt.overlap_csv
    }

    if (!defined(phenoColList) || overlap_bt.has_overlap) {
      call regenie.filter_snps as filter_snps_bt {
        input:
          bed = merge_genos.out_bed,
          bim = merge_genos.out_bim,
          fam = merge_genos.out_fam,
          samples_keep = g.right,
          prefix = group_id_bt
      }
    
      call regenie.filter_lowvar as filter_lowvar_bt {
        input:
          bed = filter_snps_bt.qc_bed
          bim = filter_snps_bt.qc_bim
          fam = filter_snps_bt.qc_fam
          pheno = split_phenotypes.bin,
          phenoColList = select_first([pheno_overlap_bt,pheno_list_bt]),
          covar = covar,
          covarColList = covarColList,
          catCovarList = catCovarList,
          bt = true,
          prefix = group_id_bt,
          R_filter_lowvar_function = R_filter_lowvar_function
      }

      call regenie.step1 as step1_bt {
        input:
          bed = filter_snps_bt.qc_bed
          bim = filter_snps_bt.qc_bim
          fam = filter_snps_bt.qc_fam
          qc_exclude = filter_lowvar_bt.lowvar_exclude,
          pheno = split_phenotypes.bin,
          phenoColList = select_first([pheno_overlap_bt,pheno_list_bt]),
          covar = covar,
          covarColList = covarColList,
          catCovarList = catCovarList,
          bt = true,
          prefix = group_id_bt
      }
    }
  }

  scatter(g in zip(cluster_traits.batches_qt, cluster_traits.samples_qt)) {
    String group_id_qt = basename(g.left, ".txt")
    String pheno_list_qt = sep(",", read_lines(g.left))

    if (defined(phenoColList)) {
      call csv.csv_overlap as overlap_qt {
        input:
          query_csv = pheno_list_qt,
          haystack_csv = select_first([phenoColList, ""])
      }
      
      String pheno_overlap_qt = overlap_qt.overlap_csv
    }

    if (!defined(phenoColList) || overlap_qt.has_overlap) {
      call regenie.filter_snps as filter_snps_qt {
        input:
          bed = merge_genos.out_bed,
          bim = merge_genos.out_bim,
          fam = merge_genos.out_fam,
          samples_keep = g.right,
          prefix = group_id_qt
     }
  
     call regenie.filter_lowvar as filter_lowvar_qt {
       input:
         bed = filter_snps_qt.qc_bed
         bim = filter_snps_qt.qc_bim
         fam = filter_snps_qt.qc_fam
         pheno = split_phenotypes.quant,
         phenoColList = select_first([pheno_overlap_qt,pheno_list_qt]),
         covar = covar,
         covarColList = covarColList,
         catCovarList = catCovarList,
         bt = false,
         prefix = group_id_qt,
         R_filter_lowvar_function = R_filter_lowvar_function
     }

     call regenie.step1 as step1_qt {
       input:
         bed = filter_snps_qt.qc_bed
         bim = filter_snps_qt.qc_bim
         fam = filter_snps_qt.qc_fam
         qc_exclude = filter_lowvar_qt.lowvar_exclude,
         pheno = split_phenotypes.quant,
         phenoColList = select_first([pheno_overlap_qt,pheno_list_qt]),
         covar = covar,
         covarColList = covarColList,
         catCovarList = catCovarList,
         bt = false,
         prefix = group_id_qt
     }
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
    File missingness_report = cluster_traits.missingness
    Array[File] diagnostic_plots = cluster_traits.diagnostics
    Array[File] pred_list = [merge_pred_list_bt.out, merge_pred_list_qt.out]
    Array[File] loco_qt = flatten(step1_qt.loco)
    Array[File] loco_bt = flatten(step1_bt.loco)
    File pheno_bt = split_phenotypes.bin
    File pheno_qt = split_phenotypes.quant
  }
}

task cluster_traits {

  input {
    File pheno
    File qc_id
    File R_batch_function
    File manifest
		String? phenoColList
    String prefix
  }

	String phenos_keep = select_first([phenoColList, ""])

  command <<<
    Rscript - <<-'CLUSTER_TRAITS'
      library(tidyverse)
      source("~{R_batch_function}")

      qc_id <- read_tsv("~{qc_id}", col_names=c("FID", "IID"))

      pheno <- read_tsv("~{pheno}") %>%
        inner_join(qc_id, by=c("eid"="IID")) %>%
        select(-FID)

      batches <- batch_by_comissingness(pheno, "~{manifest}", "~{prefix}")

      type <- map_chr(batches, ~unique(.x$type)) 

      map(c("quantitative", "binary") %>% set_names, ~which(type == .)) %>%
        iwalk(\(.x, .y) {
          paste0("batch_", .x, ".txt") %>%
            cat(file=paste0("batches_", .y, ".txt"), sep="\n")
          paste0("keep_batch_", .x, ".txt") %>%
            cat(file=paste0("keep_batches_", .y, ".txt"), sep="\n")
        })
    CLUSTER_TRAITS
  >>>

  runtime {
    memory: "64 GB"
    container: "rocker/tidyverse"
  }

  output {
    File missingness = prefix + "_cluster_missingness_summary.csv"
    Array[File] diagnostics = glob("*.png")
    Array[File] batches_qt = read_lines("batches_quantitative.txt")
    Array[File] samples_qt = read_lines("keep_batches_quantitative.txt")
    Array[File] batches_bt = read_lines("batches_binary.txt")
    Array[File] samples_bt = read_lines("keep_batches_binary.txt")
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
          select(eid, any_of(c(pheno_cols, paste0(pheno_cols, "_age"))))
      }

      split_pheno <- function(ids, out, covar) {
        pheno_split <- pheno %>%
          mutate(FID=eid) %>%
          select(FID, IID=eid, any_of(c(ids, paste0(ids, "_age"))))

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

task merge_pred_list {

  input {
    Array[File?] pred_lists
  }

  #
  # Extract the *first defined* File (safe)
  #
  Array[File] defined_files = select_all(pred_lists)

  # Ensure at least one file exists
  File first_file = defined_files[0]

  #
  # Make the merged filename:
  #   Remove "_group_<num>" from basename
  #   Append ".merged"
  #
  String base = basename(first_file)
  String cleaned = sub(base, "_group_[0-9]+", "")
  String merged_name = cleaned + ".merged"

  command <<<
    set -euo pipefail

    # Concatenate only defined files
    cat ~{sep=' ' defined_files} > ~{merged_name}
  >>>

  output {
    File out = merged_name
  }

  runtime {
    memory: "4 GB"
  }
}
