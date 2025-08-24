version 1.1

import "group_pheno.wdl" as sub

workflow regenie_step1 {

  input {
    Array[Array[File]] genos
    Array[File]+ prepared_phenotypes
    File phewas_manifest
    File? covar
    String? covarColList
    String? catCovarList
		String? phenoColList
    Float missing_thresh
  }

  scatter(p in prepared_phenotypes) {
    call sub.group_pheno {
      input:
				genos = genos,
        covar = covar,
        covarColList = covarColList,
        catCovarList = catCovarList,
        phenotype_table = p,
        phewas_manifest = phewas_manifest,
				phenoColList = phenoColList,
				missing_thresh = missing_thresh
    }
  }

  output {
    Array[File] pred_list = flatten(group_pheno.pred_list)
    Array[File] loco_qt = flatten(group_pheno.loco_qt)
    Array[File] loco_bt = flatten(group_pheno.loco_bt)
  }
}
