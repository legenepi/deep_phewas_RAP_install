version 1.0

task phenotype_preparation {

	input {
		Array[File]+ phenotype_files
		File groupings
		String phenotype_filtered_save_name
		Boolean relate_remove = true
		File? kinship_file
		Boolean IVNT = true
		String stats_save = "IVNT_RR_N_filtered_neuro_PanUKB_ancestry_stats"
		File? phewas_manifest
	}

	command <<<
		SCRIPT=`Rscript -e 'cat(system.file("extdata/scripts/association_testing","01_phenotype_preparation.R", package = "DeepPheWAS"))'` && \
		echo $SCRIPT && \
		Rscript $SCRIPT \ 
			--phenotype_filtered_save_name ~{phenotype_filtered_save_name} \
			--phenotype_files ~{sep="," phenotype_files} \
			~{"--groupings " + groupings } \
			~{true="--relate_remove" false="" relate_remove} \
			~{"--kinship_file " + kinship_file} \
			~{true="--IVNT" false="" IVNT} \
			--stats_save ~{stats_save} \
			~{"--PheWAS_manifest_overide " + phewas_manifest} \
	>>>

	output {
		File phenotype = phenotype_filtered_save_name
		File stats = stats_save
	}

	runtime {
		cpu: 1
		memory: "200 GB"
	}
}
