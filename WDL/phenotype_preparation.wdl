version development

task phenotype_preparation {

	input {
		Array[File]+ phenotype_files
		File? groupings
		String phenotype_filtered_save_name
		Boolean relate_remove
		File? kinship_file
		Boolean IVNT
		String stats_save
		File? phewas_manifest
		Int k = 120
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
		Array[File]+ out = glob("*")
	}

	runtime {
		dx_instance_type: "mem3_ssd1_v2_x32"
	}
}
