#!/usr/bin/env Rscript

DEEP_PHEWAS <- "/deep_phewas"
GENO_BASE <- "/Bulk/Imputation/UKB imputation from genotype/ukb22828_c"

'Usage:
  make_association_testing_inputs.R [--project STRING] --snp_list FILE --analysis_name STRING --out FILE

Options:
  --project STRING        RAP project ID (default: current project)
  --snp_list              SNP list file on the RAP 
  --analysis_name STRING  RAP output analysis name
  --out FILE              JSON file to create
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(tidyverse))
suppressMessages(library(jsonlite))
source("make_inputs_functions.R")

args <- docopt(doc)

phenotype_files <- c("composite_phenotypes.RDS", "concepts.RDS", "all_dates.RDS",
                     "range_ID.RDS", "phecodes.RDS", "PQP.RDS", "formula_phenotypes.RDS",
                     "data_field_phenotypes.RDS") %>%
  map_chr(~paste(DEEP_PHEWAS, ., sep="/") %>% get_file_id(args$project))

options_files <- list("association_testing.phenotype_preparation.kinship_file" = "related_callrate",
                      "association_testing.phenotype_preparation.groupings" = "inputs/ancestry_panUKB",
                      "association_testing.covariates" = "inputs/covariates_plink")

plink_phenotype_files <- c("AFR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz",
                           "AMR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz",
                           "CSA_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz",
                           "EAS_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz",
                           "EUR_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz",
                           "MID_IVNT_RR_N_filtered_neuro_PanUKB_ancestry.gz") %>%
  map_chr(~paste(DEEP_PHEWAS, "phenotypes", ., sep="/") %>% get_file_id(args$project))

snps_bgis <- 
  
inputs <- list("association_testing.phenotype_preparation.phenotype_files" = phenotype_files,
               "association_testing.phenotype_preparation.relate_remove" = TRUE,
               "association_testing.prepare_phenotypes" = FALSE,
               "association_testing.snp_list" = get_file_id(args$snp_list, args$project))
               "association_testing.phenotype_preparation.stats_save" =
                 "IVNT_RR_N_filtered_neuro_PanUKB_ancestry_stats",
               "association_testing.analysis_name" = args$analysis_name,
               "association_testing.phenotype_preparation.IVNT" = TRUE,
               "association_testing.phenotype_preparation.phenotype_filtered_save_name" =
                 "IVNT_RR_N_filtered_neuro_PanUKB_ancestry",
               "association_testing.PLINK_association_testing.phenotypes" = plink_phenotype_files,
               "association_testing.extracting_snps.bgis": "Array[File]+ (optional, default = [\"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c10_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c6_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c17_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c9_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c14_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c21_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c4_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cXY_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c11_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c20_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c19_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c7_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c12_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c22_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c18_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c5_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c15_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c8_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c2_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c16_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cX_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c3_b0_v3.bgen.bgi\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c13_b0_v3.bgen.bgi\"])",
               
               
               "association_testing.extracting_snps.bgens": "Array[File]+ (optional, default = [\"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c10_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c6_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c17_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c9_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c14_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c21_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c4_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cXY_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c11_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c20_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c19_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c7_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c12_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c22_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c18_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c5_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c15_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c8_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c2_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c16_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_cX_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c3_b0_v3.bgen\", \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c13_b0_v3.bgen\"])",
               "association_testing.phewas_manifest": "File? (optional, default = \"dx://project-GJbvyPjJy3Gy01jz4x8bXzgv:/deep_phewas/inputs/PheWAS_manifest.csv\")"

write_json(inputs, args$out, pretty=TRUE, auto_unbox=TRUE)