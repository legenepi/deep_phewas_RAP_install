#!/usr/bin/env Rscript

'Usage:
  make_inputs_phenotype_generation --data_files=FILES --GPC=FILE --GPP=FILE --hesin_diag=FILE --HESIN=FILE --hesin_oper=FILE --death_cause=FILE --death=FILE --out=FILE [ --save_location=STRING --exclusions=FILE --king_coef=FILE --PheWAS_manifest_overide=FILE --concept_codes=FILE --PQP_codes=FILE --composite_phenotype_map_overide=FILE ]

Options:
    --data_files=FILES                        Comma separated full file paths of the data files that will be formatted and concatenated.
    --GPC=FILE                                Full path of the primary care clinical data from UK Biobank.
    --GPP=FILE                                Full path of the primary care prescription data from UK Biobank.
    --hesin_diag=FILE                         Full file path of hesin_diag file from UK-Biobank.
    --HESIN=FILE                              Full file path of HESIN file from UK-Biobank.
    --hesin_oper=FILE                         Full file path of the hesin_oper file from UK-Biobank.
    --death_cause=FILE                        Full file path of the death_cause file from UK-Biobank.
    --death=FILE                              Full file path of the death file from UK-Biobank contains date of death info.
    --save_location=STRING                    Full file path for the common folder to save created files
    --exclusions=FILE                         File containing individuals to be excluded from the analysis
    --king_coef=FILE                          Full file path of related data file containing King coefficient scores for related ID pairs
    --PheWAS_manifest_overide=FILE            Alternative PheWAS_manifest file
    --concept_codes=FILE                      Zip file containing concept code lists
    --PQP_codes=FILE                          Zip file containing alternative primary care code lists
    --composite_phenotype_map_overide=FILE    Composite_phenotype_map file
    --out=FILE                                JSON input file to create
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(tidyverse))
suppressMessages(library(jsonlite))

PROJECT_ID <- Sys.getenv("PROJECT_ID")
INPUTS <- Sys.getenv("INPUTS")

source("R/make_inputs_functions.R")

args <- docopt(doc)

file_inputs <- list(phenotype_generation.data_preparation.death=args$death,
                    phenotype_generation.data_preparation.death_cause=args$death_cause,
                    phenotype_generation.data_preparation.GPC=args$GPC,
                    phenotype_generation.data_preparation.GPP=args$GPP,
                    phenotype_generation.data_preparation.HESIN=args$HESIN,
                    phenotype_generation.data_preparation.hesin_diag=args$hesin_diag,
                    phenotype_generation.data_preparation.hesin_oper=args$hesin_oper,
                    phenotype_generation.data_preparation.king_coef=args$king_coef) %>%
    map(~get_file_id(., PROJECT_ID))

optional_inputs <- list(phenotype_generation.minimum_data.exclusions=args$exclusions,
                        phenotype_generation.concept_codes=args$concept_codes,
                        phenotype_generation.PQP_codes=args$PQP_codes,
                        phenotype_generation.phewas_manifest=args$PheWAS_manifest_overide,
                        phenotype_generation.composite_phenotype_map_overide=args$composite_phenotype_map_overide) %>%
    discard(is.null) %>%
    map(~get_upload_id(., PROJECT_ID, INPUTS))

inputs <- c(list(phenotype_generation.minimum_data.files=map(args$data_files, ~get_file_id(., PROJECT_ID))),
            file_inputs,
            optional_inputs)

write_json(inputs, args$out, pretty=TRUE, auto_unbox=TRUE)
