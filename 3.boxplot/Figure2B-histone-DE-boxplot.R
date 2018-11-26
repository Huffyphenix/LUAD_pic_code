
# gene list ---------------------------------------------------------------
.libPaths("E:/library")
genelist <- c("CBX2","EZH2")
genelist <- c("TP53")
genelist <- c("E2F1","E2F2","E2F3","E2F5","SOX4","NME2","TP63","TP53","CBX7")
genelist <- c("CCNA2","CCNB1",
              "CCNB2","CCNE1",
              "CCNE2","CDK1",
              "CDKN2A","CDC25C",
              "MCM2","MCM4","MCM6","MCM7")
genelist <- c("JARID2","AEBP2","EED","SET","SUZ12","RBBP4","RBBP7","EZH1") #PRC2 complex
genelist <- c("RYBP","RING1","RNF2","BMI1","PCGF2","PCGF1","PHC1","PHC2","PHC3") #PRC1 complex
genelist <- c("CBX2","CBX4","CBX6","CBX7","CBX8") # other CBX

de_path <- "G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
# de_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"

# data_path_3 <- "S:/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
data_path_3 <- "G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
library(magrittr)
TF_DE_info <- read.table(file.path(de_path,"NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T) %>%
  dplyr::rename("Gene_id"="gene_id")
progene_DE_info <- read.table(file.path(de_path,"NOISeq_DE_ProGene_FC2_cpm_30"),sep = '\t',header = T) 
rbind(TF_DE_info,progene_DE_info) -> all_DE_info
# all gene no filter: cpm>1
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil
# load exp ----------------------------------------------------------------

exp <- readr::read_tsv("G:/data/TCGA/lung_DE/CBX24678/all_genes_exp.confirm")
exp <- readr::read_tsv(file.path(de_path,"NOISeq_DE_ProGene_FC2_cpm_30.exp.xls"))

# filter ------------------------------------------------------------------

library(magrittr)
exp %>%
  dplyr::filter(gene_id %in% genelist) %>%
  tidyr::gather(-gene_id,key="sample",value="Expression") %>%
  # .[-1,] %>%
  # dplyr::rename("PRDM5_RSEM"="value") %>%
  dplyr::mutate(Group=substr(sample,6,7)) %>%
  dplyr::mutate(Group=ifelse(Group=="01","Tumor","Normal")) %>%
  dplyr::mutate(log2Exp=log2(Expression)) -> genelist_exp



all_gene_nofil %>%
  dplyr::filter(gene_id %in% genelist) %>%
  dplyr::select(gene_id,log2FC) %>%
  dplyr::mutate(title=paste(gene_id," (log2FC = ",round(log2FC,2),")",sep=""))-> genelist_FC

genelist_exp %>%
  dplyr::group_by(gene_id) %>%
  dplyr::do(broom::tidy(t.test(Expression ~Group,data=.))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::select(gene_id,fdr) %>%
  dplyr::mutate(label=paste0("FDR = ",signif(fdr, 3)))-> genelist_stage_pvalue
# draw pic ----------------------------------------------------------------

# use ggpubr --------------------------------------------------------------
library(ggplot2)
genelist_exp %>%
  dplyr::inner_join(genelist_FC,by="gene_id") %>%
  dplyr::arrange(Group) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(ID=substr(sample,1,4)) -> df
b <- runif(nrow(df), -0.2, 0.2)

comp_list <- list(c("1","2"))
df %>%
  ggpubr::ggboxplot(x = "Group", y = "log2Exp",
                    color = "Group", palette = "npg" #add = "jitter",
                    ) +
  facet_wrap(~ title, strip.position = "bottom", nrow = 1,scales = "free") +
  geom_point(aes(x=as.numeric(Group)+b,y=log2Exp,color=Group)) +
  geom_line(aes(x=as.numeric(Group)+b,y=log2Exp,group=ID),linetype="11",color="grey") +
  scale_color_manual(
    values = c("#1E90FF","#EE6363")
  )+
  ylab("log2(mRNA Exp)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white",colour = "white"),
        strip.text = element_text(size = 12)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif") +
  scale_x_discrete(breaks = c(1,2),
                     labels = c("Normal","Tumor"),
                     expand = c(0.2,0.2)) -> p;p


ggsave("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure2B.DE_histone_boxplot.pdf",device = "pdf",width = 5,height = 3)
ggsave("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure2B.DE_histone_boxplot.tiff",device = "tiff",width = 5,height = 3)

ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary/Figure S1.TP53_boxplot.pdf",device = "pdf",width = 4,height = 5)
ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary/Figure S1.cell_cycle_genes_boxplot.pdf",device = "pdf",width = 8,height = 5)
ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure1/Figure1B.DE_histone_boxplot.pdf",device = "pdf",width = 6,height = 5)

ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure S2.PRC2_boxplot.pdf",device = "pdf",width = 6,height = 5)
ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure S2.PRC2_boxplot.tiff",device = "tiff",width = 6,height = 5)

ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure S2.PRC1_boxplot.pdf",device = "pdf",width = 6,height = 5)
ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure S2.PRC1_boxplot.tiff",device = "tiff",width = 6,height = 5)

ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure3/Figure S3.CBX_boxplot.pdf",device = "pdf",width = 5,height = 4)

ggsave("S:/study/ENCODE-TCGA-LUAD/CBX4678/Figure2B.DE_histone_boxplot.pdf",device = "pdf",width = 10,height = 4)
ggsave("S:/study/ENCODE-TCGA-LUAD/CBX4678/Figure2B.DE_histone_boxplot.tiff",device = "tiff",width = 10,height = 4)

# by ggplot2 --------------------------------------------------------------

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
  ) -> p;p

ann_text <- data.frame(Group = "Tumor",Expression = 2000,lab = genelist_stage_pvalue$label,
                       gene_id = factor(genelist_stage_pvalue$gene_id,levels = genelist_stage_pvalue$gene_id))
p + geom_text(data = ann_text,aes(label=ann_text$lab))
ggsave("F:/?ҵļ?????/ENCODE-TCGA-LUAD/Figure/Figure2/Figure2B.DE_histone_boxplot.pdf",device = "pdf",width = 10,height = 4)
