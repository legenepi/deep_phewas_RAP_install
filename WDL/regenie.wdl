version 1.1

task filter_genos {

  input {
    File bed
    File bim
    File fam
    File pheno
  }

  String out = basename(bed, ".bed")

  command <<<
    zcat "~{pheno}" | awk 'NR > 1 { print $1, $1}' > samples.keep
    plink2 --bed "~{bed}" --bim "~{bim}" --fam "~{fam}" --keep samples.keep --make-bed --out "~{out}_filt"
  >>>

  output {
    File out_bed = "~{out}_filt.bed"
    File out_bim = "~{out}_filt.bim"
    File out_fam = "~{out}_filt.fam"
  }

  runtime {
    memory: "12 GB"
    container: "us.gcr.io/broad-dsde-methods/plink2-alpha"
  }
}

task merge_genos {

  input {
    Array[File] beds
    Array[File] bims
    Array[File] fams
    String prefix
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
    memory: "100 GB"
    container: "us.gcr.io/broad-dsde-methods/plink2-alpha"
  }
}

task filter_snps {

  input {
    File bed
    File bim
    File fam
    File? samples_keep
    String prefix
  }

  String out = prefix + "_qc_pass"

  command <<<
    plink2 --bed "~{bed}" --bim "~{bim}" --fam "~{fam}" \
      ~{"--keep " + samples_keep} \
      --geno 0.1 \
      --hwe 1e-15 \
      --mac 100 \
      --maf 0.01 \
      --mind 0.1 \
      --no-id-header \
      --out "~{out}" \
      --write-samples \
      --write-snplist
  >>>

  output {
    File qc_id = out + ".id"
    File qc_snplist = out + ".snplist"
  }

  runtime {
    memory: "64 GB"
    container: "us.gcr.io/broad-dsde-methods/plink2-alpha"
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

  String out = prefix + "_step1"

  command <<<
    BED="~{bed}"
    regenie\ 
      --step 1 \
      --bed "${BED/.bed/}" \
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
    memory: "64 GB"
    disks: "1000 GB"
    container: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
  }
}
