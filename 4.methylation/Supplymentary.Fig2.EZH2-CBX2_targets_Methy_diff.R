
# package prepair ---------------------------------------------------------

library(magrittr,ggplot2)

# data path ---------------------------------------------------------------

methy_data_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/"
genelist_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"

# data manage -------------------------------------------------------------

genelist_pro <- readr::read_tsv(file.path(genelist_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
genelist_TF <- readr::read_tsv(file.path(genelist_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))
methy <- readr::read_rds(file.path(methy_data_path,"pan33_allgene_methy_diff.simplification.rds.gz"))
methy_cor <- readr::read_rds(file.path(methy_data_path,"pancan34_all_gene_exp-cor-meth.rds.gz"))

# gene list data ----------------------------------------------------------

methy %>%
  dplyr::filter(cancer_types=='LUAD') %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) -> LUAD_gene_methy

methy_cor %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) -> LUAD_gene_methy_cor

LUAD_gene_methy %>%
  dplyr::filter(diff >=0) -> LUAD_gene_methy_diff_0

LUAD_gene_methy_cor %>%
  dplyr::filter(symbol %in% LUAD_gene_methy_diff_0$symbol) %>%
  dplyr::filter(spm <= -0.4) -> genes_significantly_cor_and_diff.0

LUAD_gene_methy %>%
  dplyr::filter(symbol %in% genes_significantly_cor_and_diff.0$symbol)
