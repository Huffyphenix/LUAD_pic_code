library(magrittr)
data_path<-"S:/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"S:/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"
immune_histone<-read.table(file.path(data_path,"immune_histone.txt"),sep = "\t",header = T)
hist(immune_histone$EZH2_karyon)

#####################################
#correlation analysis
#pearson 
broom::tidy(
    cor.test(immune_histone$CBX2_cytoplsm,immune_histone$EZH2_cytoplsm,method = "pearson")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high) ->CBX2_EZH2_cytoplsm
rownames(CBX2_EZH2_cytoplsm)="CBX2_EZH2_cytoplsm"

broom::tidy(
  cor.test(immune_histone$CBX2_cytoplsm,immune_histone$EZH2_karyon,method = "pearson")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::select(estimate,p.value,fdr,conf.low,conf.high) ->CBX2_EZH2_karyon
rownames(CBX2_EZH2_karyon)="CBX2_EZH2_karyon"

rbind(CBX2_EZH2_cytoplsm,CBX2_EZH2_karyon) ->cbx2_ezh2_result
cbx2_ezh2_result %>%
  write.table(file=file.path(result_path,"pearson-result"),quote=F,sep="\t")

#########################heatmap
rownames(immune_histone) <-immune_histone$sample_type %>% paste(c(1:nrow(immune_histone)),sep=".")

immune_histone %>%
  dplyr::select(CBX2_cytoplsm,EZH2_karyon,EZH2_cytoplsm) %>%
  pheatmap::pheatmap(cluster_cols = F)

##########################stage
oneway.test(EZH2_karyon~Stage,data = immune_histone)
boxplot(EZH2_karyon~Stage,data = immune_histone)

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

ggsave(filename = "immune_histone-stage.pdf", plot = p, device = "pdf", path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2", 
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
  dplyr::mutate(sample_type=ifelse(sample_type=="TA","Normal(TA)",sample_type)) %>%
  dplyr::arrange(sample_type) %>%
  ggpubr::ggboxplot(x = "sample_type", y = "PositiveRateXStainingIntensity",
                    color = "sample_type", pallete = "npg",add = "jitter",
                    facet.by = "group") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 2.3)->p;p
ggsave(filename = "immune_histone.pdf", plot = p, device = "pdf", path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",
       width = 8, height = 5)


save.image(file = file.path(result_path,"immune_histone.rdata"))
