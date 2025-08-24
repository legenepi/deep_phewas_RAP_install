version 1.1

import "WDL/group_pheno.wdl"

workflow regenie_step2_sv {

  input {
		Array[Array[File]] genos
    File phewas_manifest
    File snp_list
    String analysis_name
    File? covariates
    File? phenotype_inclusion_file
    Array[File]+ phenotypes
    String? covarColList
    String? catCovarList
    File pred_list_bt
    Array[File] loco_bt
    File pred_list_qt
    Array[File] loco_qt
  }

	Array[Pair[File, Array[File]]] crossed = cross(phenotypes, genos)

  scatter (pg in crossed) {
    call group_pheno.split_phenotypes as split_pheno {
      input:
        phenotype_table = pg.left,
        phewas_manifest = phewas_manifest,
        covar = covariates
    }
         
    call regenie_association_testing as regenie_bt {
      input:
        bgen = pg.right[0],
        bgi = pg.right[1],
        sample = pg.right[2],
				extract = snp_list,
        pheno = split_pheno.bin,
        covar = covariates,
        covarColList = covarColList,
        catCovarList = catCovarList,
        pred_list = pred_list_bt,
        loco = loco_bt,
        bt = true
    }

    call regenie_association_testing as regenie_qt {
      input:
        bgen = pg.right[0],
        bgi = pg.right[1],
        sample = pg.right[2],
				extract = snp_list,
        pheno = split_pheno.quant,
        covar = covariates,
        covarColList = covarColList,
        catCovarList = catCovarList,
        pred_list = pred_list_qt,
        loco = loco_qt,
        bt = false
    }
  }
	
	scatter (p in phenotypes) {
    call merge_sv {
			input:
				prefix = analysis_name,
			  pheno = p,
				results = flatten([ regenie_bt.results, regenie_qt.results ]),
				dict = regenie_bt.dict
		}		
	}

  output {
    Array[File] sv_results = merge_sv.merged
  }
}

task regenie_association_testing {

  input {
    File bgen
    File bgi
    File sample
    File pheno
    File? covar
    File? exclude
    File? extract
    String? covarColList
    String? catCovarList
    File pred_list
    Array[File] loco
    Boolean bt
  }

  String base = basename(bgen, '.bgen') 
  String out = "base" + "~{if bt then '_bin' else '_quant'}"

  command <<<
    ln -s ~{sep=" " loco} .
    regenie \
      --step 2 \
      --bgen "~{bgen}" \
			--bgi "~{bgi}" \
			--sample "~{sample}" \
      ~{"--exclude " + exclude } \
      ~{"--extract " + extract } \
      ~{"--covarFile " + covar} \
      ~{"--covarColList " + covarColList} \
      ~{"--catCovarList " + catCovarList} \
      --phenoFile "~{pheno}" \
      --pred "~{pred_list}" \
      ~{true="--bt" false="--qt" bt} \
      --bsize 400 \
      --minMAC 3 \
      --threads=16 \
      --nauto 23 \
      --no-split \
      --out "~{out}"
  >>>

  output {
    File results = "~{out}.regenie"
    File dict = "~{out}.regenie.Ydict"
  }

  runtime {
    dx_instance_type: "mem3_ssd1_v2_x16"
    container: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
  }
}

task merge_sv {

  input {
    String prefix
    String pheno
    Array[File] results
    Array[File] dict
  }

  String out = "~{prefix}_sv_~{pheno}.regenie.gz"

  command <<<
    awk 'NR == FNR {
      pheno[$1]=$2
      npheno++
      next
    }
    NR == npheno + 1 {
      for (i in pheno)
        gsub(i, pheno[i], $0)
      print
    }
    FNR > 2 { print }' "~{dict[0]}" `echo "~{sep=' ' results}" | grep -o "[^ ]*~{pheno}[^ ]*"` | gzip > ~{out}
  >>>
  
  output {
    File merged = out
  }

  runtime {
    cpu: 2
    memory: "16 GB"
  }
}
