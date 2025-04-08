version 1.1

import "phenotype_preparation.wdl"
import "group_pheno.wdl" as sub

workflow association_stage1 {

  input {
    File? covar
    String? covarColList
    String? catCovarList
    Array[File]+ phenotype_files
    File? groupings
    String phenotype_filtered_save_name
    Boolean relate_remove
    File? kinship_file
    Boolean IVNT
    String stats_save
    File phewas_manifest
  }

  call phenotype_preparation.phenotype_preparation {
    input:
      phenotype_files = phenotype_files,
      groupings = groupings,
      phenotype_filtered_save_name = phenotype_filtered_save_name,
      relate_remove = relate_remove,
			kinship_file = kinship_file,
      IVNT = IVNT,
      stats_save = stats_save,
      phewas_manifest = phewas_manifest
  }

  scatter(p in phenotype_preparation.out) {
    call sub.group_pheno {
      input:
        covar = covar,
        covarColList = covarColList,
        catCovarList = catCovarList,
        phenotype_table = p,
        phewas_manifest = phewas_manifest
    }
  }

  output {
    Array[File] pred_list = flatten(group_pheno.pred_list)
    Array[File] loco_qt = flatten(group_pheno.loco_qt)
    Array[File] loco_bt = flatten(group_pheno.loco_bt)
  }
}
