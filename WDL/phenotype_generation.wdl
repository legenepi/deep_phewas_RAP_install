version 1.0

import "data_wrangling.wdl"
import "phenotype_creation.wdl"

workflow phenotype_generation {

	input {
		Array[File]+ tab_data
		File GPC
		File GPP
		File hesin_diag
		File HESIN
		File hesin_oper
		File death
		File death_cause
		File? king_coef
		File? exclusions
		String save_loc = "."
		File? phewas_manifest
		File? concept_codes
		File? PQP_codes
		File? composite_phenotype_map_overide
	}

	call data_wrangling.minimum_data {
			input:
				files = tab_data,
				exclusions = exclusions,
				save_loc = save_loc
	}
	
	call data_wrangling.data_preparation {
		input:
			min_data = minimum_data.out,
			GPC = GPC,
			GPP = GPP,
			hesin_diag = hesin_diag,
			HESIN = HESIN,
			hesin_oper = hesin_oper,
			death = death,
			death_cause = death_cause,
			king_coef = king_coef,
			save_loc = save_loc
	}

	call phenotype_creation.phecode_generation {
		input:
			health_data = data_preparation.health_records,
			sex_info = data_preparation.combined_sex,
			control_exclusions = data_preparation.control_exclusions
	}

	call phenotype_creation.data_field_phenotypes {
		input:
			min_data = minimum_data.out,
			phewas_manifest = phewas_manifest
	}

	call phenotype_creation.creating_concepts {
		input:
			GPP = data_preparation.GP_P,
			health_data = data_preparation.health_records,
			phewas_manifest = phewas_manifest,
			code_list = concept_codes
	}

	call phenotype_creation.primary_care_quantitative_phenotypes {
		input:
			GPC = data_preparation.GP_C,
			DOB = data_preparation.DOB,
			phewas_manifest = phewas_manifest,
			code_list = PQP_codes
	}

	call phenotype_creation.formula_phenotypes {
		input:
			min_data = minimum_data.out,
			data_field_phenotypes = data_field_phenotypes.out,
			sex_info = data_preparation.combined_sex,
			phewas_manifest = phewas_manifest
	}

	call phenotype_creation.composite_phenotypes {
		input:
			phecode_file = phecode_generation.phecode_file,
			range_ID_file = phecode_generation.range_ID_file,
			concept_file = creating_concepts.concept_file,
			all_dates_file = creating_concepts.all_dates_file,
			data_field_file = data_field_phenotypes.out,
			PQP_file = primary_care_quantitative_phenotypes.out,
			formula_file = formula_phenotypes.out,
			control_populations = data_preparation.control_populations,
			composite_phenotype_map_overide = composite_phenotype_map_overide
	}

	output {
		File min_data = minimum_data.out
		Array[File]+ out = data_preparation.out
		Array[File]+ phecode = phecode_generation.out
		File data_field_pheno = data_field_phenotypes.out
		Array[File]+ concept = creating_concepts.out
		File pqp = primary_care_quantitative_phenotypes.out
		File formula = formula_phenotypes.out
		File composite = composite_phenotypes.out
	}
}
