library(magrittr)
data_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"
immune_histone<-read.table(file.path(data_path,"immune_histone.txt"),sep = "\t",header = T) %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="T","1Tumor","2Normal"))
hist(immune_histone$EZH2_karyon)

#####################################
#correlation analysis
#pearson 
immune_histone %>%
  dplyr::filter(sample_type=="1Tumor") -> immune_histone.T
immune_histone %>%
  dplyr::filter(sample_type=="2Normal") -> immune_histone.N
broom::tidy(
    cor.test(immune_histone.T$CBX2_cytoplsm,immune_histone.T$EZH2_cytoplsm,method = "pearson")) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high) %>%
  dplyr::mutate(x=1,y=1.8,sample_type="1Tumor") %>%
  dplyr::mutate(label=paste("Pearson Cor = ",round(estimate,2),"\n P.value = ",format(p.value,scientific=TRUE,digit=2),"     \n",
                            "n = ",nrow(immune_histone.T),
                            "                    ",sep="")) ->CBX2_EZH2_cytoplsm.T
broom::tidy(
  cor.test(immune_histone.N$CBX2_cytoplsm,immune_histone.N$EZH2_cytoplsm,method = "pearson")) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high) %>%
  dplyr::mutate(x=0.41,y=0.8,sample_type="2Normal") %>%
  dplyr::mutate(label=paste("Pearson Cor = ",round(estimate,2),"\n P.value = ",format(p.value,scientific=TRUE,digit=2),"    \n",
                            "n = ",nrow(immune_histone.N),
                            "                     ",sep=""))->CBX2_EZH2_cytoplsm.N
rbind(CBX2_EZH2_cytoplsm.T,CBX2_EZH2_cytoplsm.N) ->CBX2_EZH2_cytoplsm
facet_names <- list(
  '1Tumor'="Tumor",
  '2Normal'="Normal"
)
facet_labeller <- function(variable,value){
  return(facet_names[value])
}
immune_histone %>%
  ggplot(aes(x=EZH2_cytoplsm,y=CBX2_cytoplsm)) +
  geom_point(aes(color = sample_type)) +
  geom_smooth(se = FALSE, fullrange=TRUE, color = "#039BE5") +
  facet_wrap(~sample_type,scales = "free",labeller=facet_labeller) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none"
  ) +
  geom_text(data=CBX2_EZH2_cytoplsm,aes(x=x,y=y,label=label),hjust=0.5) +
  labs(
    x = "EZH2_cytoplsm",
    y = "CBX2_cytoplsm"
  ) -> p1;p1
ggsave(filename = "EZH2-cytoplsm_CBX2-cytoplsm_immunehistochemistry_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 4,height = 3)
ggsave(filename = "EZH2-cytoplsm_CBX2-cytoplsm_immunehistochemistry_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 4,height = 3)

broom::tidy(
  cor.test(immune_histone.T$CBX2_cytoplsm,immune_histone.T$EZH2_karyon,method = "pearson")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high)  %>%
  dplyr::mutate(x=1,y=1.8,sample_type="1Tumor") %>%
  dplyr::mutate(label=paste("Pearson Cor = ",round(estimate,2),"\n P.value = ",format(p.value,scientific=TRUE,digit=2),"     \n",
                            "n = ",nrow(immune_histone.T),
                            "                     ",sep="")) ->CBX2_EZH2_karyon.T
broom::tidy(
  cor.test(immune_histone.N$CBX2_cytoplsm,immune_histone.N$EZH2_karyon,method = "pearson")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high)  %>%
  dplyr::mutate(x=0.41,y=0.8,sample_type="2Normal") %>%
  dplyr::mutate(label=paste("Pearson Cor = ",round(estimate,2),"\n P.value = ",format(p.value,scientific=TRUE,digit=2),"   \n",
                            "n = ",nrow(immune_histone.N),
                            "                    ",sep="")) ->CBX2_EZH2_karyon.N

rbind(CBX2_EZH2_karyon.T,CBX2_EZH2_karyon.N) ->CBX2_EZH2_karyon

immune_histone %>%
  ggplot(aes(x=EZH2_karyon,y=CBX2_cytoplsm)) +
  geom_point(aes(color = sample_type)) +
  geom_smooth(se = FALSE, fullrange=TRUE, color = "#039BE5") +
  facet_wrap(~sample_type,scales = "free",labeller = facet_labeller) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none"
  ) +
  geom_text(data=CBX2_EZH2_karyon,aes(x=x,y=y,label=label),hjust=0.5) +
  labs(
    x = "EZH2_karyon",
    y = "CBX2_cytoplsm"
  ) -> p2;p2
ggsave(filename = "EZH2-karyo_CBX2-cytoplsm_immunehistochemistry_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 4,height = 3)
ggsave(filename = "EZH2-karyo_CBX2-cytoplsm_immunehistochemistry_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 4,height = 3)


library(grid)
library(gridExtra)
library(scales)
p2 + theme(axis.title.y = element_blank()) ->p2
pdf(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.pdf"),width = 8,height = 3)
tiff(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.tiff"),width = 8,height = 3,units ="in",res=400)

grid.arrange(p1, p2, ncol=2,nrow=1,widths=c(1,1), heights=c(1))
dev.off()



#########################heatmap
rownames(immune_histone) <-immune_histone$sample_type %>% paste(c(1:nrow(immune_histone)),sep=".")

immune_histone %>%
  dplyr::select(CBX2_cytoplsm,EZH2_karyon,EZH2_cytoplsm) %>%
  pheatmap::pheatmap(cluster_cols = F)

##########################stage
oneway.test(EZH2_karyon~stage,data = immune_histone)
boxplot(EZH2_karyon~stage,data = immune_histone)

library(ggplot2)
comp_list <- list(c("Normal(TA)", "stage II"), c("stage II", "stage III"))
immune_histone %>%
  dplyr::mutate(stage=as.character(stage)) %>%
  dplyr::mutate(stage=ifelse(is.na(stage),"Normal(TA)",stage)) %>%
  dplyr::arrange(stage) %>%
  tidyr::gather(-c("sample","sample_type","stage"),key="group",value="PositiveRateXStainingIntensity") %>%
  # dplyr::filter(group=="EZH2_karyon") %>%
  ggpubr::ggboxplot(x = "stage", y = "PositiveRateXStainingIntensity",
                    color = "stage", palette = "npg", add = "jitter",
                    facet.by = "group") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 3) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(2, 2.5))->p;p

ggsave(filename = "immune_histone-stage.pdf", plot = p, device = "pdf", path = "F:/?ҵļ?????/ENCODE-TCGA-LUAD/Figure/Figure2", 
       width = 8, height = 5)

# immune_histone %>%
#   dplyr::filter(sample_type=="T") %>%
#   ggpubr::ggboxplot(x = "stage", y = "CBX2_cytoplsm",  color = c(2,3), pallete = "jco"  ) +
#   ggpubr::stat_compare_means(method = "anova")->p;p
# ggsave(filename = "CBX2_cytoplsm_stage.pdf", plot = p, device = "pdf", path = result_path, width = 5, height = 5)
# 
# immune_histone %>%
#   dplyr::filter(sample_type=="T") %>%
#   ggpubr::ggboxplot(x = "stage", y = "EZH2_cytoplsm",  color = c(2,3), pallete = "jco"  ) +
#   ggpubr::stat_compare_means(method = "anova")->p;p
# ggsave(filename = "EZH2_cytoplsm_stage.pdf", plot = p, device = "pdf", path = result_path, width = 5, height = 5)

##########################Tumor and Tumor adjcent
immune_histone %>%
  tidyr::gather(-c("sample","sample_type","stage"),key="group",value="PositiveRateXStainingIntensity") %>%
  dplyr::mutate(sample_type=as.character(sample_type)) %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="TA","2","1")) %>%
  dplyr::arrange(sample_type) ->df
df %>%
  ggpubr::ggboxplot(x = "sample_type", y = "PositiveRateXStainingIntensity",
                    color = "sample_type", pallete = "npg",add = "jitter",
                    facet.by = "group") +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("2","1"),
                   labels=c("Normal","Tumor")) +
  theme(
    axis.title.x = element_blank()
  ) +
  ylab("PositiveRate X tainingIntensity") +
  ggpubr::stat_compare_means(label.y = 2.3,method = "kruskal.test")->p;p
ggsave(filename = "immune_histone.tiff", plot = p, device = "tiff", path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",
       width = 8, height = 3)


save.image(file = file.path(result_path,"immune_histone.rdata"))
