
# package prepair ---------------------------------------------------------

library(magrittr,ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
# data path ---------------------------------------------------------------

methy_data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/"
genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
chip_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets"
data_path<- "H:/data"

# output data -------------------------------------------------------------
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# data manage -------------------------------------------------------------
genelist <- readr::read_tsv(file.path(chip_path,"common-targets-180426-new","all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(Class == "Down") %>%
  .$gene_id.x
genelist <- readr::read_tsv(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea","gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description %in% "PPAR signaling pathway")%>%
  .$SYMBOL
cell_cycle_relate <- c("Cell cycle","Oocyte meiosis","DNA replication",
                       "Homologous recombination","p53 signaling pathway",
                       "Progesterone-mediated oocyte maturation")
genelist <- readr::read_tsv(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea","gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description %in% cell_cycle_relate)%>%
  .$SYMBOL

Animal_TF <-  readr::read_tsv(file.path(data_path,"AnimalTFDB","Homo_sapiens_transcription_factors_gene_list.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls")) %>%
  .$Symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Animal_TF %>%
  dplyr::filter(! Entrez_ID %in% enzyme_lsit$ENTREZID) -> Animal_TF
Animal_TF$Entrez_ID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db) -> Animal_TF_entrez

# genelist_pro <- readr::read_tsv(file.path(genelist_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
# genelist_TF <- readr::read_tsv(file.path(genelist_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))
methy <- readr::read_rds(file.path(methy_data_path,"pan33_allgene_methy_diff.simplification.rds.gz"))
methy_cor <- readr::read_rds(file.path(methy_data_path,"pancan34_all_gene_exp-cor-meth.rds.gz"))

tcga_geneid <- readr::read_tsv("F:/我的坚果云/ENCODE-TCGA-LUAD/TCGA_gene_info/TCGA_all_gene_id.txt") %>%
  dplyr::rename("symbol"="gene_id")

# gene list data ----------------------------------------------------------

methy %>%
  dplyr::filter(cancer_types=='LUAD') %>%
  tidyr::unnest() %>%
  # dplyr::inner_join(tcga_geneid,by="symbol") %>%
  # dplyr::filter(symbol %in% genelist$SYMBOL) %>%
  dplyr::filter(symbol %in% genelist) -> LUAD_gene_methy

methy_cor %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() %>%
  dplyr::inner_join(tcga_geneid,by="symbol") %>%
  # dplyr::filter(entrez_id %in% genelist$ENTREZID) %>%
  dplyr::filter(symbol %in% genelist) -> LUAD_gene_methy_cor

LUAD_gene_methy %>%
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  readr::write_tsv(file.path(out_path_sup,"correlation-table/cellcycle_pathway_methy_diff_cor.tsv"))

### For downregulate genes ------
LUAD_gene_methy %>%
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  dplyr::filter(spm<=-0.3 & diff>0) %>%
  readr::write_tsv(file.path(out_path_fig,"Figure4/Figure5","genes_regulate_by_methy.tsv"))

LUAD_gene_methy %>%
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  dplyr::filter(spm<=-0.3 & diff>0) -> LUAD_gene_methy.sig_gene

# draw pic ----------------------------------------------------------------
LUAD_gene_methy_cor %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::arrange(spm) %>% .$symbol -> cor_rank.genesymbol
LUAD_gene_methy_cor %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::rename("value"="spm") %>%
  dplyr::mutate(group="Cor.") -> LUAD_gene_methy_cor.pic
LUAD_gene_methy %>%
  dplyr::rename("value"="diff","logfdr" = "fdr") %>%
  dplyr::mutate(group="Diff. (T - N)") %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::select(-direction) -> LUAD_gene_methy.pic

library(ggplot2)
library(grid)
CPCOLS <- c("red", "white", "#1C86EE")

LUAD_gene_methy_cor.pic %>%
  dplyr::select(-entrez_id) %>%
  rbind(LUAD_gene_methy.pic) %>%
  dplyr::mutate(value=signif(value,3)) %>%
  ggplot(aes(x=group,y=symbol)) +
  geom_tile(aes(fill = value),color="white") +
  scale_y_discrete(limit = cor_rank.genesymbol) +
  guides(size = guide_legend(title.position = "left",
                             title.theme = element_text(angle = 90))) +
  scale_fill_gradient2(
    name = "Diff./Cor.", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
  ) +
  geom_text(aes(label=value)) +
  ylab("Symbol") +
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
ggsave(file.path(out_path_fig,"Figure1","PPAR_methy_Cor-diff-gsca.pdf"),device = "pdf",width = 4,height = 4)
ggsave(file.path(out_path_fig,"Figure1","PPAR_methy_Cor-diff-gsca.tiff"),device = "tiff",width = 4,height = 4)

ggsave(file.path(out_path_fig,"Figure4/Figure5","Figure5B.methy_Cor-diff-gsca.pdf"),device = "pdf",width = 4,height = 6)
ggsave(file.path(out_path_fig,"Figure4/Figure5","Figure5B.methy_Cor-diff-gsca.tiff"),device = "tiff",width = 4,height = 6)

### For upregulate genes -------
LUAD_gene_methy %>%
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  dplyr::filter(spm<=-0.3 & diff<0) %>%
  readr::write_tsv(file.path(out_path_fig,"Figure4/Figure5","genes_regulate_by_methy.tsv"))

LUAD_gene_methy %>%
  dplyr::inner_join(LUAD_gene_methy_cor,by="symbol") %>%
  dplyr::filter(spm<=-0.3 & diff<0) -> LUAD_gene_methy.sig_gene

LUAD_gene_methy_cor %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::arrange(spm) %>% .$symbol -> cor_rank.genesymbol

LUAD_gene_methy_cor %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::rename("value"="spm") %>%
  dplyr::mutate(group="Cor.") -> LUAD_gene_methy_cor.pic
LUAD_gene_methy %>%
  dplyr::rename("value"="diff","logfdr" = "fdr") %>%
  dplyr::mutate(group="Diff. (T - N)") %>%
  dplyr::filter(symbol %in% LUAD_gene_methy.sig_gene$symbol) %>%
  dplyr::select(-direction) -> LUAD_gene_methy.pic
CPCOLS <- c("#00688B", "#00BFFF")

LUAD_gene_methy_cor.pic %>%
  dplyr::select(-entrez_id) %>%
  rbind(LUAD_gene_methy.pic) %>%
  dplyr::mutate(value=signif(value,3)) %>%
  ggplot(aes(x=group,y=symbol)) +
  geom_tile(aes(fill = value),color="white") +
  scale_y_discrete(limit = cor_rank.genesymbol) +
  guides(size = guide_legend(title.position = "left",
                             title.theme = element_text(angle = 90))) +
  scale_fill_gradient2(
    name = "Diff./Cor.", #"Methylation diff (T - N)",
    low = CPCOLS[2],
    high = CPCOLS[1],
    breaks = c(-0.6,-0.4,-0.2,0)
  ) +
  geom_text(aes(label=value)) +
  ylab("Symbol") +
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


ggsave(file.path(out_path_fig,"Figure1","cellcycle.methy_Cor-diff-gsca.pdf"),device = "pdf",width = 4,height = 4)
ggsave(file.path(out_path_fig,"Figure1","cellcycle.methy_Cor-diff-gsca.tiff"),device = "tiff",width = 4,height = 4)

