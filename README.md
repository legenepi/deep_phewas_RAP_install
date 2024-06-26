# Overview
These scripts install and run [DeepPheWAS](https://github.com/Richard-Packer/DeepPheWAS) on the [UK Biobank Research Analysis Platform (RAP)](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform).

## Requirements

* R
* [`dx-toolkit`](https://github.com/dnanexus/dx-toolkit)
* Java
* You must be logged into the RAP (`dx login`)

## Configuration

You should only have to edit the `options.config` file before running the installation and association testing scripts.

## Installing DeepPheWAS in your RAP project

1. In `options.config` set your RAP project ID (PROJECT_ID), the install location (PROJECT_DIR) and the filename of your most up-to-date RAP dataset (DATASET).
2. Run `./install_workflows.sh` which will:
    - Run a RAP job that creates a docker image in your project containing the DeepPheWAS package, plink2 and dependencies.
    - Download [DNANexus dxCompiler](https://documentation.dnanexus.com/developer/building-and-executing-portable-containers-for-bioinformatics-software/dxcompiler) required to compile WDL workflows to DNANexus workflows to run on the RAP.
    - Compile and install the [WDL](https://github.com/openwdl/wdl) DeepPheWAS workflows in your RAP project.

## Extracting the UK Biobank data fields that provide the input data for phenotype generation

Run `./extract_fields.sh` which will launch several `table-exporter` jobs to extract:

* Participant data fields
* Hospital episode statistics (HES)
* Death registry data
* Primary care data

## Phenotype generation

1. In `options.config` set:
    - Any sample exclusions (exclusions=exclusions.txt)
    - Optionally, specify any files for construction of additional bespoke phenotypes (phewas_manifest, concept_codes, PQP_codes, composite_phenotype_map_overide).

You can leave these options blank and the default set of phenotypes will be generated (see [DeepPheWAS documentation](https://richard-packer.github.io/DeepPheWAS_site/) for details of options and creation of additional phenotypes).

2. Run `./make_inputs_phenotype_generation.sh` which will make the .json configuration file specifying the inputs required for the phenotype generation step.
3. Run `./run_phenotype_generation.sh` which will submit the RAP job to generate the phenotypes.

## Phenotype preparation

1. In `options.config` set:
    - Whether to remove related samples (related_remove=true/false)
    - Whether to rank inverse-normal transform phenotypes (IVNT=true/false)
    - Any groupings of samples, e.g. by ancestry (groupings=ancestry_panUKB)
2. Run `./make_inputs_phenotype_preparation.sh` which will make the .json configuration file specifying the inputs required for the phenotype preparation step.
3. Run `./run_phenotype_preparation.sh` which will submit the RAP job to create the phenotype tables to be used in association testing.

## Association testing

1. Prepare a file containing the list of SNPs to be tested with the format below:
```
chromosome,rsid,group_name,coded_allele,non_coded_allele,graph_save_name
1,rs12097169,LYSMD1,A,C,rs12097169_A_C
1,rs4845556,TDRKH,A,G,rs4845556_A_G
2,rs59985551,EFEMP1,T,C,rs59985551_T_C
```
3. In `options.config` set:
    - The SNP list file created above (snp_list=mysnps.csv)
    - The output name you want for the analysis (analysis_name=myproject)
4. Run `./run_association_testing.sh` which will create the .json file specifying inputs for the job and submits the job to the RAP.

The outputs should be found in your DeepPheWAS installation directory under `results/<analysis_name>`, which you would then download with `dx download`.
