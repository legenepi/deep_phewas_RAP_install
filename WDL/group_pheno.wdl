version 1.1

import "regenie_step1.wdl"

workflow group_pheno {

  input { 
    Array[Array[File]] genos = [
    ["dx://file-FxXZzV0JkF65g2vX9Vx8jkkZ", "dx://file-FxXZxjQJkF636yj326QqXPV5", "dx://file-Gbq1G7jJ8Jj31pbkFYyXPy4G"],
    ["dx://file-FxXb030JkF6GbKYb9VP89bpf", "dx://file-FxXZxp8JkF67GV2z7v5ZqFzZ", "dx://file-Gbq18jQJ8Jj3Gzy9jZ8fbGvG"],
    ["dx://file-FxXb0X0JkF6JV6g0Pp0gvZq8", "dx://file-FxXZxpQJkF69vjv312xj39xz", "dx://file-Gbq189QJ8JjPGyvJ6Gy04k4Y"],
    ["dx://file-FxXb0y8JkF67pkJx13BV0J7G", "dx://file-FxXZxq0JkF65bZBj13QpBZFZ", "dx://file-Gbq18J0J8Jj06j03ggygjBx9"],
    ["dx://file-FxXb1F0JkF68b43G0pf356x6", "dx://file-FxXZxq8JkF61JXY78Y2xp2Z6", "dx://file-Gbq18V0J8Jj05QVFkQK1Kj5f"],
    ["dx://file-FxXb1bjJkF6B6F1Z9VYv7fxv", "dx://file-FxXZxqjJkF6PpqqK0p2gZxvb", "dx://file-Gbq181QJ8JjJK73pkf3FzzYP"],
    ["dx://file-FxXb238JkF6GX0PQ9VB1Gqjv", "dx://file-FxXZxv8JkF6457pxPkfgJ60Y", "dx://file-Gbq1860J8JjFb6yFVJ2g2Zg8"],
    ["dx://file-FxXb2JjJkF6BxKfyJYkyFFKQ", "dx://file-FxXZxvQJkF67Z48f0pjq53Q8", "dx://file-Gbq17XjJ8Jj5p478zZ3kgJQ1"],
    ["dx://file-FxXb2j0JkF6B5yVBPkVvV376", "dx://file-FxXZxx0JkF6BXzY51BBbz1f9", "dx://file-Gbq17f8J8JjPyq9XPjKqbBpK"],
    ["dx://file-FxXb30QJkF641JG196p7gj6B", "dx://file-FxXZxZjJkF64ZZKv1B2V667F", "dx://file-Gbq1Jb0J8Jj9kbBQGGb4BY4v"],
    ["dx://file-FxXb390JkF67KVkZPpF2PYz9", "dx://file-FxXZxb0JkF69X2FQ19v6jXVQ", "dx://file-Gbq1Gj8J8Jj8kPb6XX0bK9Zk"],
    ["dx://file-FxXb3b0JkF60ZXbP9Vkf47b2", "dx://file-FxXZxbQJkF6G2F900kkp0YY4", "dx://file-Gbq1GvQJ8Jj5XpjXkk1K1YJG"],
    ["dx://file-FxXb3y0JkF62xzjX7vv52ZyG", "dx://file-FxXZxbjJkF6GbX1619P5YBj9", "dx://file-Gbq1J0jJ8Jj8kPb6XX0bK9b3"],
    ["dx://file-FxXb45jJkF64ZZKv1B2V675y", "dx://file-FxXZxf8JkF6BXzY51BBbz1bG", "dx://file-Gbq1J4jJ8JjPYxv26zvXGBbz"],
    ["dx://file-FxXb4K0JkF68b43G0pf3574X", "dx://file-FxXZxfjJkF67Z48f0pjq53PG", "dx://file-Gbq1GG0J8JjPqpF49z2jY9xb"],
    ["dx://file-FxXb4f8JkF6JV6g0Pp0gvb8Q", "dx://file-FxXZxg0JkF68b43G0pf356Y0", "dx://file-Gbq1GQ0J8Jj82xg0jP8pzPbj"],
    ["dx://file-FxXb4x0JkF6B6F1Z9VYv7gG6", "dx://file-FxXZxgQJkF66qV4X0p7KBgJG", "dx://file-Gbq1GZQJ8JjB75b4x6VPpqYQ"],
    ["dx://file-FxXb54QJkF6JV6g0Pp0gvb9p", "dx://file-FxXZxgjJkF6B6gJg1B6qXGjf", "dx://file-Gbq1FzjJ8Jj3k8QpJjP19XgX"],
    ["dx://file-FxXb5BjJkF6GbX1619P5YFVb", "dx://file-FxXZxj8JkF6B0kfB44f4PpVK", "dx://file-Gbq1G3QJ8Jj0gXbGf770VZzP"],
    ["dx://file-FxXb5VjJkF6BxKfyJYkyFFzB", "dx://file-FxXZxk0JkF6GbKYb9VP89bYB", "dx://file-Gbq18v0J8Jj7QjG0FV80y2K0"],
    ["dx://file-FxXb5fjJkF6K07FF9X4Q2PG0", "dx://file-FxXZxk8JkF67GV2z7v5ZqFzQ", "dx://file-Gbq1900J8Jj045YZ81Vz5ZXG"],
    ["dx://file-FxXb5pQJkF6PpqqK0p2gZzv1", "dx://file-FxXZxkjJkF641JG196p7ggbF", "dx://file-Gbq18ZQJ8JjGjQ9gQ8gQ3g98"]
  ]
    File? covar
    String? covarColList
    String? catCovarList
    File phenotype_table
    File phewas_manifest
    Int k_high = 120
    Int k_low = 10
    Float missing_thresh = 0.15
  }

  String prefix = basename(phenotype_table, ".tsv.gz")

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
      pheno = split_phenotypes.quant,
      covar = covar,
      qc_id = filter_snps.qc_id,
      qc_snplist = filter_snps.qc_snplist,
      covarColList = covarColList,
      catCovarList = catCovarList,
      bt = false,
      prefix = basename(split_phenotypes.quant, ".txt")
  }

  call calc_dist {
    input:
      pheno = split_phenotypes.bin,
      qc_id = filter_snps.qc_id
  }
  
  call cluster_traits {
    input:
      pheno = split_phenotypes.bin,
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
        pheno = split_phenotypes.bin,
        phenoColList = g.left,
        covar = covar,
        qc_id = g.right,
        qc_snplist = filter_snps.qc_snplist,
        covarColList = covarColList,
        catCovarList = catCovarList,
        bt = true,
        prefix = basename(split_phenotypes.bin, ".txt") + "_" + basename(g.right, ".ids")
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

  String out_bin = basename(phenotype_table, ".tsv.gz") + "_bin.txt"
  String out_quant = sub(out_bin, "bin", "quant")

  command <<<
    Rscript - <<-SPLIT_PHENO
      library(tidyverse)
      manifest <- read_csv("~{phewas_manifest}",
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

      split_pheno(id_bin, "~{out_bin}")
      split_pheno(id_quant, "~{out_quant}") 
    SPLIT_PHENO
  >>>

  output {
    File bin = out_bin
    File quant = out_quant
  }

  runtime {
    memory: "200 GB"
    container: "rocker/tidyverse"
  }
}

task calc_dist {

  input {
    File pheno
    File qc_id
  }

  String dist_out = basename(pheno, ".txt") + "_dist.RDS"

  command <<<
    Rscript - <<-CALC_DIST
      install.packages("parallelDist")
      library(parallelDist)
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
    memory: "200 GB"
    container: "rocker/tidyverse"
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

  String groups_out = basename(pheno, ".txt") + "_group"
  Int n_groups = k_high + k_low
  Array[Int] group_idx = range(n_groups)

  command <<<
    Rscript - <<-CLUSTER_TRAITS
      install.packages("usedist")
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
                      select(FID, IID, all_of(.x %>% pull(pheno))) %>%
                      filter(if_any(.x %>% pull(pheno), ~!is.na(.))))
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
              paste0("~{groups_out}_", ., ".ids"), col_names = FALSE))
      walk(c(~{k_high + 1}:~{n_groups})-1, ~write_tsv(pheno_low_missing.list[[.]] %>% select(FID, IID),
              paste0("~{groups_out}_", ., ".ids"), col_names = FALSE))
    CLUSTER_TRAITS
  >>>

  runtime {
    memory: "64 GB"
    container: "rocker/tidyverse"
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

  String pred_list_merged = sub(basename(select_first(pred_lists)), "_group*_", "")

  command <<<
    cat "~{sep=' ' pred_lists}" > ~{pred_list_merged}
  >>>

  runtime {
    memory: "4 GB"
  }

  output {
    File out = pred_list_merged
  }
}
