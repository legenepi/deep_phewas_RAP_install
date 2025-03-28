version development

import "regenie_step1.wdl"

workflow group_pheno {

	input { 
		Array[Array[File]] genos
		File? covar
		String? covarColList
		String? catCovarList
		File phenotype_table
		File phewas_manifest
		Int k_high = 120
		Int k_low = 10
		Float missing_thresh = 0.15
	}

	String prefix = basename(phenotype_table, ".txt")

	scatter(geno in genos) {
		call regenie_step1.filter_genos {
			input:
				bed = geno[0],
				bim = geno[1],
				fam = geno[2],
				pheno = phenotype_table
		}
	}

	call regenie_step1.merge_genos {
		input:
			beds = filter_genos.out_bed,
			bims = filter_genos.out_bim,
			fams = filter_genos.out_fam,
			prefix = prefix
	}

	call regenie_step1.filter_snps {
		input:
			bed = merge_genos.out_bed,
			bim = merge_genos.out_bim,
			fam = merge_genos.out_fam,
			prefix = prefix
	}

	call split_phenotypes {
		input:
			phenotype_table = phenotype_table,
			phewas_manifest = phewas_manifest
	}

	call regenie_step1.step1 as step1_qt {
		input:
			bed = merge_genos.out_bed,
			bim = merge_genos.out_bim,
			fam = merge_genos.out_fam,
			pheno = split_phenotypes.out_quant,
			covar = covar,
			qc_id = filter_snps.qc_id,
			qc_snplist = filter_snps.qc_snplist,
			covarColList = covarColList,
      catCovarList = catCovarList,
      bt = false,
			prefix = basename(split_phenotypes.out_quant, ".txt")
	}

	call calc_dist {
		input:
			pheno = split_phenotypes.out_bin,
			qc_id = filter_snps.qc_id
	}
	
	call cluster_traits {
		input:
			pheno = split_phenotypes.out_bin,
			dist = calc_dist.out,
			k_high = k_high,
			k_low = k_low,
			missing_thresh = missing_thresh,
			qc_id = filter_snps.qc_id
	}
	
	scatter(g in zip(cluster_traits.pheno_groups, cluster_traits.group_ids)) {
		call regenie_step1.step1 as step1_bt {
			input:
				bed = merge_genos.out_bed,
				bim = merge_genos.out_bim,
				fam = merge_genos.out_fam,
				pheno = split_phenotypes.out_bin,
				phenoColList = g.left,
				covar = covar,
				qc_id = g.right,
				qc_snplist = filter_snps.qc_snplist,
				covarColList = covarColList,
				catCovarList = catCovarList,
				bt = true,
				prefix = basename(split_phenotypes.out_bin, ".txt") + "_" + basename(g.right, ".ids")
		}
	}

	call merge_pred_list {
		input:
			pred_lists = step1_bt.pred_list
	}

	output {
		Array[File] pred_list = [step1_qt.pred_list, merge_pred_list.out]
		Array[File] loco_qt = step1_qt.loco
		Array[File] loco_bt = flatten(step1_bt.loco)
	}
}

task split_phenotypes {

	input {
		File phenotype_table
		File phewas_manifest
	}

	command <<<
		Rscript <<-SPLIT_PHENO
			library(tidyverse)
			manifest <- read_csv("PheWAS_manifest.csv",
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

			pheno <- read_tsv("~{phenotype_table}", n_max = 1000)

			split_pheno <- function(ids, out) {
				pheno %>%
					mutate(FID=eid) %>%
					select(FID, IID=eid, any_of(ids)) %>%
					write_tsv(out)
			}

			split_pheno(id_bin, str_replace("~{phenotype_table}", ".tsv.gz", "_bin.txt"))
			split_pheno(id_quant, str_replace("~{phenotype_table}", ".tsv.gz", "_quant.txt")) 
		SPLIT_PHENO
	>>>

	output {
		File out_bin = glob("*_bin.txt")
		File out_quant = glob("*_quant.txt")
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x32"
		container: "rocker/r-ver:4.3.1"
	}
}

task calc_dist {

	input {
		File pheno
		File qc_id
	}

	String dist_out = sub(pheno, ".txt", "_dist.RDS")

	command <<<
		Rscript - <<-CALC_DIST
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
		dx_instance_type: "mem3_ssd1_v2_x32"
		container: "rocker/r-ver:4.3.1"
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
		Int k_high
		Int k_low
		Float missing_thresh
	}

	String groups_out = sub(pheno, ".txt", "_group")
	Int n_groups = k_high + k_low
	Array[Int] group_idx = range(n_groups)

	command <<<
		Rscript - <<-CLUSTER_TRAITS
			library(tidyverse)
			
			qc_id <- read_tsv("~{qc_id}", col_names=c("FID", "IID"))

			pheno_bin <- read_tsv("~{pheno}") %>%
				inner_join(qc_id, by=c("FID", "IID"))

			pheno_bin.missing <- pheno_bin %>%
				summarise(across(starts_with("P"), ~sum(is.na(.)/length(.)))) %>%
				pivot_longer(everything(), names_to = "pheno", values_to = "missing")

			cluster_phenos <- function(x.dist, k=120, method="ward.D2") {
				hclust(x.dist, method=method) %>%
					cutree(k=k) %>%
					enframe(name="pheno", value="group") %>%
					mutate(group=factor(group)) %>%
					inner_join(pheno_bin.missing, by="pheno") %>%
					group_by(group) %>%
					group_map(~pheno_bin %>%
											select(c(FID, IID, .x$pheno)) %>%
											filter(if_any(.x$pheno, ~!is.na(.))))
			}

			calc_missing <- function(x) {
				summarise(x, across(c(-FID, -IID), ~sum(is.na(.)/length(.)))) %>%
					pivot_longer(everything(), names_to = "pheno", values_to = "missing")
			}

			pheno_dist <- readRDS("~{dist}")

			pheno_high_missing <- usedist::dist_subset(pheno_dist, pheno_bin.missing %>%
				filter(missing >= ~{missing_thresh}) %>%
				pull(pheno))

			pheno_high_missing.list <- pheno_high_missing %>%
				cluster_phenos(k=~{k_high})

			pheno_high_missing_summary <- pheno_high_missing.list %>%
				map_df(calc_missing, .id = "group")

			pheno_low_missing <- usedist::dist_subset(pheno_dist, pheno_bin.missing %>%
				filter(missing < ~{missing_thresh}) %>%
				pull(pheno))

			pheno_low_missing.list <- pheno_low_missing %>%
				cluster_phenos(k=~{k_low})

			pheno_low_missing_summary <- pheno_low_missing.list %>%
				map_df(calc_missing, .id = "group")

			pheno_low_missing_groups <- pheno_low_missing_summary %>%
				mutate(group = sample(c(~{k_high + 1}:~{n_groups}), n(), replace = TRUE) %>% factor)

			pheno_groups <- bind_rows(list(pheno_high_missing_summary, pheno_low_missing_groups))

			pheno_groups_out <- pheno_groups %>%
				select(group, pheno) %>%
				group_by(group) %>%
				summarise(pheno=paste(pheno, collapse = ",")) %>%
				ungroup %>%
				mutate(group=as.integer(group)) %>%
				arrange(group) %>%
				pull(pheno)

			cat(pheno_groups_out, file="~{groups_out}.txt", sep="\n")

			walk(c(1:~{k_high})-1, ~write_tsv(pheno_high_missing.list[[.]] %>% select(FID, IID),
							paste0("~{groups_out}_, ., ".ids"), col_names = FALSE))
			walk(c(~{k_high + 1}:~{n_groups})-1, ~write_tsv(pheno_low_missing.list[[.]] %>% select(FID, IID),
							paste0("~{groups_out}_, ., ".ids"), col_names = FALSE))
		CLUSTER_TRAITS
	>>>

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x8"
		container: "rocker/r-ver:4.3.1"
	}

	output {
		Array[String] pheno_groups = read_lines(groups_out)
		Array[File] group_ids = suffix(".ids", prefix("group_", group_idx))
	}
}

task merge_pred_list {

	input {
		Array[File] pred_lists
	}

	String pred_list_merged = sub(select_first(pred_lists), "_group*_", "")

	command <<<
		cat "~{sep=' ' pred_lists}" > ~{pred_list_merged}
	>>>

	runtime {
		    dx_instance_type: "mem1_ssd1_v2_x2" 
	}

	output {
		File out = pred_list_merged
	}
}
