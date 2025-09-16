version 1.1

struct Step1 {
  File pheno
  File pred_list
  Array[File] loco
}

workflow regenie_step2_collapsing {

  input {
    Array[Array[File]] genos
    File? gene_list
    String analysis_name
    File? covariates
    String? phenoColList
    Array[File] phenotypes_bt
    Array[File] phenotypes_qt
    String? covarColList
    String? catCovarList
    Array[File] pred_list_bt
    Array[Array[File]] loco_bt
    Array[File] pred_list_qt
    Array[Array[File]] loco_qt
    File? exclude
    File annot
    File setlist
    File mask
    Array[Float] aaf_bins
    String joint_tests
    Float maxaff
    String tests
  }

  Array[Pair[File, Array[File]]] pred_loco_bt = zip(pred_list_bt, loco_bt)
  Array[Pair[File, Pair[File, Array[File]]]] pheno_pred_loco_bt = zip(phenotypes_bt, pred_loco_bt)

  scatter(ppl in pheno_pred_loco_bt) {
    Step1 step1_bt = Step1 { pheno: ppl.left, pred_list: ppl.right.left, loco: ppl.right.right }
  }

  Array[Pair[Step1, Array[File]]] crossed_bt = cross(step1_bt, genos)

  scatter (pg in crossed_bt) {
    call regenie_association_testing as regenie_bt {
      input:
        bgen = pg.right[0],
        bgi = pg.right[1],
        sample = pg.right[2],
        exclude = exclude,
        pheno = pg.left.pheno,
        phenoColList = phenoColList,
        covar = covariates,
        covarColList = covarColList,
        catCovarList = catCovarList,
        pred_list = pg.left.pred_list,
        loco = pg.left.loco,
        bt = true,
        prefix = analysis_name,
        annot = annot,
        setlist = setlist,
        mask = mask,
        aaf_bins = aaf_bins,
        joint_tests = joint_tests,
        maxaff = maxaff,
        tests = tests,
        gene_list = gene_list
    }
  }

  Array[Pair[File, Array[File]]] pred_loco_qt = zip(pred_list_qt, loco_qt)
  Array[Pair[File, Pair[File, Array[File]]]] pheno_pred_loco_qt = zip(phenotypes_qt, pred_loco_qt)

  scatter(ppl in pheno_pred_loco_qt) {
    Step1 step1_qt = Step1 { pheno: ppl.left, pred_list: ppl.right.left, loco: ppl.right.right }
  }

  Array[Pair[Step1, Array[File]]] crossed_qt = cross(step1_qt, genos)

  scatter (pg in crossed_qt) {
    call regenie_association_testing as regenie_qt {
      input:
        bgen = pg.right[0],
        bgi = pg.right[1],
        sample = pg.right[2],
        exclude = exclude,
        pheno = pg.left.pheno,
        phenoColList = phenoColList,
        covar = covariates,
        covarColList = covarColList,
        catCovarList = catCovarList,
        pred_list = pg.left.pred_list,
        loco = pg.left.loco,
        bt = false,
        prefix = analysis_name,
        annot = annot,
        setlist = setlist,
        mask = mask,
        aaf_bins = aaf_bins,
        joint_tests = joint_tests,
        maxaff = maxaff,
        tests = tests,
        gene_list = gene_list
    }
  }
  
  output {
    Array[File] sv_results = flatten([regenie_bt.results, regenie_qt.results])
    Array[File] sv_dict = flatten([regenie_bt.dict, regenie_qt.dict])
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
    String? covarColList
    String? catCovarList
    String? phenoColList
    File pred_list
    Array[File] loco
    Boolean bt
    File annot
    File setlist
    File mask
    Array[Float] aaf_bins
    String joint_tests
    Float maxaff
    String tests
    String prefix
    File? gene_list
    File? extract_sets
  }

  String pheno_base = basename(pheno, '.txt')
  String bgen_base = basename(bgen, '.bgen') 
  String out = "~{prefix}_~{pheno_base}_~{bgen_base}"

  command <<<
    ln -s ~{sep=" " loco} .
    regenie \
      --step 2 \
      --bgen "~{bgen}" \
      --bgi "~{bgi}" \
      --sample "~{sample}" \
      ~{"--exclude " + exclude } \
      ~{"--covarFile " + covar} \
      ~{"--covarColList " + covarColList} \
      ~{"--catCovarList " + catCovarList} \
      --phenoFile "~{pheno}" \
      ~{"--phenoColList " + phenoColList} \
      ~{"--extract-sets " + gene_list } \
      --pred "~{pred_list}" \
      ~{true="--bt" false="--qt" bt} \
      --bsize 200 \
      --anno-file "~{annot}" \
      --set-list "~{setlist}" \
      --mask-def "~{mask}" \
      --nauto 23 \
      --aaf-bins "~{sep=',' aaf_bins}" \
      --joint "~{joint_tests}" \
      --vc-maxAAF "~{maxaff}" \
      --vc-tests "~{tests}" \
      --rgc-gene-p \
      --threads=16 \
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
