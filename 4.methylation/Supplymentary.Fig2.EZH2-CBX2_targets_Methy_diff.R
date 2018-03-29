
# package prepair ---------------------------------------------------------

library(magrittr,ggplot2)

# data path ---------------------------------------------------------------

methy_data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/"
genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"

# output data -------------------------------------------------------------
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

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
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  readr::write_tsv(file.path(out_path_sup,"CBX2_EZH2_targets_methy_diff_cor.tsv"))


LUAD_gene_methy_cor %>%
  dplyr::arrange(spm) %>% .$symbol -> cor_rank.genesymbol
# draw pic ----------------------------------------------------------------
LUAD_gene_methy_cor %>%
  dplyr::rename("value"="spm") %>%
  dplyr::mutate(group="Spearman Cor") -> LUAD_gene_methy_cor.pic
LUAD_gene_methy %>%
  dplyr::rename("value"="diff","logfdr" = "fdr") %>%
  dplyr::mutate(group="Methylation diff (T - N)") %>%
  dplyr::select(-direction) -> LUAD_gene_methy.pic

CPCOLS <- c("red", "white", "blue")
LUAD_gene_methy_cor.pic %>%
  rbind(LUAD_gene_methy.pic) %>%
  ggplot(aes(x=group,y=symbol)) +
  geom_point(aes(size=logfdr,color=value)) +
  scale_y_discrete(limit = cor_rank.genesymbol) +
  scale_color_gradient2(
    name = "Diff/Cor", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1],
    breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
  ) +
  scale_size_continuous(
    name = "-log10(Pvalue)"
  ) +
  theme(#legend.position = "bottom",
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = "white", colour = "black") ,
        plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(out_path_fig,"Figure4","Figure4B.methy_Cor-diff-gsca.pdf"),device = "pdf",width = 6,height = 9)


