version development

task filter_genos {

	input {
		File bed
		File bim
		File fam
		File pheno
	}

	String out = basename(bed, ".bed")

	command <<<
		plink2 --bed "~{bed}" --bim "~{bim}" --fam "~{fam}" --keep "~{pheno}" --make-bed --out "~{out}_filt"
	>>>

	output {
		File out_bed = "~{out}_filt.bed"
		File out_bim = "~{out}_filt.bim"
		File out_fam = "~{out}_filt.fam"
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x8"
		container: "pgscatalog/plink2"
	}
}

task merge_genos {

	input {
		Array[File] beds
		Array[File] bims
		Array[File] fams
		String prefix
	}
	
	parameter_meta {
		beds : "stream"
	}

	command <<<
		cat "~{write_lines(beds)}" | sed -e 's/.bed//g' > merge_list.txt
		plink2 --pmerge-list merge_list.txt bfile --make-bed --out "~{prefix}"
	>>>

	output {
		File out_bed = prefix + ".bed"
		File out_bim = prefix + ".bim"
		File out_fam = prefix + ".fam"
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x16"
		container: "pgscatalog/plink2"
	}
}

task filter_snps {

	input {
		File bed
		File bim
		File fam
		String prefix
	}

	command <<<
		plink2 --bed "~{bed}" --bim "~{bim}" --fam "~{fam}" \
			--geno 0.1 \
			--hwe 1e-15 \
			--mac 100 \
			--maf 0.01 \
			--mind 0.1 \
			--no-id-header \
			--out qc_pass \
			--write-samples \
			--write-snplist
	>>>

	output {
		File qc_id = prefix + "_qc_pass.id"
		File qc_snplist = prefix + "_qc_pass.snplist"
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x8"
		container: "pgscatalog/plink2"
	}
}

task step1 {

	input {
		File bed
		File bim
		File fam
		File pheno
		File? covar
		File? qc_id
		File? qc_snplist
		String? covarColList
		String? catCovarList
		String? phenoColList
		Boolean bt
		String prefix
	}

	String bed_path = sub(bed, ".bed", "")
	String out = prefix + "_step1"

	parameter_meta {
		bed : "stream"
	}

	command <<<
		regenie\ 
			--step 1 \
			--bed "~{bed_path}" \
			~{"--extract " + qc_snplist} \
			~{"--keep " + qc_id} \
			~{"--covarFile " + covar} \
			~{"--covarColList " + covarColList} \
			~{"--catCovarList " + catCovarList} \
			--phenoFile "~{pheno}" \
			~{"--phenoColList " + phenoColList} \
			--bsize 1000 \
			--lowmem \
			--lowmem-prefix . \
			~{true="--bt" false="--qt" bt} \
			--use-relative-path \
			--out "~{out}" \
			--threads 8
	>>>

	output {
		File pred_list = "~{out}_pred.list"
		Array[File] loco = glob("~{out}_*.loco")
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x8"
		container: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
	}
}
