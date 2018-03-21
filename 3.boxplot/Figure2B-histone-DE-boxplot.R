
# gene list ---------------------------------------------------------------

genelist <- c("CBX2","EZH2","UHRF1")

de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/�������data"
TF_DE_info <- read.table(file.path(de_path,"FC2","NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T) %>%
  dplyr::rename("Gene_id"="gene_id")
progene_DE_info <- read.table(file.path(de_path,"FC2","NOISeq_DE_ProGene_FC2_cpm_30"),sep = '\t',header = T) 
rbind(TF_DE_info,progene_DE_info) -> all_DE_info
# load exp ----------------------------------------------------------------

exp <-readr::read_tsv("H:/data/TCGA/lung_DE/CBX24678/all_genes_exp.confirm")


# filter ------------------------------------------------------------------

library(magrittr)
exp %>%
  dplyr::filter(gene_id %in% genelist) %>%
  tidyr::gather(-gene_id,key="sample",value="Expression") %>%
  # .[-1,] %>%
  # dplyr::rename("PRDM5_RSEM"="value") %>%
  dplyr::mutate(Group=substr(sample,6,7)) %>%
  dplyr::mutate(Group=ifelse(Group=="01","Tumor","Normal")) -> genelist_exp

all_DE_info %>%
  dplyr::filter(Gene_id %in% genelist) %>%
  dplyr::select(Gene_id,log2FC) -> genelist_FC

genelist_exp %>%
  dplyr::group_by(gene_id) %>%
  dplyr::do(broom::tidy(t.test(Expression ~Group,data=.))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::select(gene_id,fdr) %>%
  dplyr::mutate(label=paste0("FDR = ",signif(fdr, 3)))-> genelist_stage_pvalue
# draw pic ----------------------------------------------------------------

library(ggplot2)
genelist_exp$gene_id %>%
  plyr::revalue(c(CBX2="CBX2 (Log2FC = 2.55)",
                  EZH2="EZH2 (Log2FC = 2.75)",
                  UHRF1="UHRF1 (Log2FC = 3.8)")) ->genelist_exp$gene_id
genelist_exp%>%
  ggplot2::ggplot(mapping=aes(x=Group,y=Expression,color=Group)) +
  geom_boxplot() +
  geom_point(aes(x=Group,y=Expression,color=Group), position = "jitter") +
  facet_grid(~gene_id) +
  theme(
    axis.line = element_line(color = "black"),
    panel.background  = element_rect(fill = "white", color = "grey"),
    panel.grid = element_line(colour = "grey"),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20)
  ) -> p

ann_text <- data.frame(Group = "Tumor",Expression = 2000,lab = genelist_stage_pvalue$label,
                       gene_id = factor(genelist_stage_pvalue$gene_id,levels = genelist_stage_pvalue$gene_id))
p + geom_text(data = ann_text,aes(label=ann_text$lab))
ggsave("F:/�ҵļ����/ENCODE-TCGA-LUAD/Figure/Figure2/Figure2B.DE_histone_boxplot.pdf",device = "pdf",width = 10,height = 4)