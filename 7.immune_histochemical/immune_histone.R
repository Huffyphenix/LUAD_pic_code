.libPaths("E:/library")
library(magrittr)
# HOME -----
data_path<-"Z:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"

# E Zhou -----
data_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"

# HUST ----
data_path<-"G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"

immune_histone<-read.table(file.path(data_path,"immune_histone.txt"),sep = "\t",header = T) %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="T","Tumor","Normal")) %>%
  dplyr::mutate(EZH2_cytoplsm=ifelse(is.na(EZH2_cytoplsm),0,EZH2_cytoplsm)) %>%
  dplyr::mutate(EZH2_karyon=ifelse(is.na(EZH2_karyon),0,EZH2_karyon)) %>%
  dplyr::mutate(CBX2_cytoplsm=ifelse(is.na(CBX2_cytoplsm),0,CBX2_cytoplsm)) %>%
  dplyr::mutate(CBX2_karyon=ifelse(is.na(CBX2_karyon),0,CBX2_karyon)) %>%
  dplyr::mutate(EZH2=EZH2_cytoplsm+EZH2_karyon) %>%
  dplyr::mutate(CBX2=CBX2_cytoplsm+CBX2_karyon)
hist(immune_histone$EZH2_karyon)

# Preleminary test to check the test assumptions
shapiro.test(immune_histone$EZH2_cytoplsm) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(immune_histone$EZH2_karyon) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(immune_histone$CBX2_cytoplsm) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(immune_histone$CBX2_karyon) # p-value < 0.05, don't follow a normal distribution.

#####################################
#correlation analysis
#pearson 
immune_histone %>%
  dplyr::filter(sample_type=="Tumor") -> immune_histone.T
immune_histone %>%
  dplyr::filter(sample_type=="Normal") -> immune_histone.N

human_read <- function(.x){
  if (.x > 0.1) {
    .x %>% signif(digits = 2) %>% toString()
  } else if (.x < 0.1 && .x > 0.001 ) {
    .x %>% signif(digits = 1) %>% toString()
  } else {
    .x %>% format(digits = 2, scientific = TRUE)
  }
}

broom::tidy(
    cor.test(immune_histone.T$CBX2,immune_histone.T$EZH2,method = "kendall")) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=1,y=1.8,sample_type="1Tumor",n=nrow(immune_histone.T), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2.T
broom::tidy(
  cor.test(immune_histone.N$CBX2_cytoplsm,immune_histone.N$EZH2_cytoplsm,method = "kendall")) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=0.41,y=0.8,sample_type="2Normal",n=nrow(immune_histone.N), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2_cytoplsm.N
rbind(CBX2_EZH2_cytoplsm.T,CBX2_EZH2_cytoplsm.N) %>%
  dplyr::as.tbl() ->CBX2_EZH2_cytoplsm
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
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  scale_color_manual(
    values = c("#EE6363","#00C5CD"),
    labels = CBX2_EZH2_cytoplsm$label
  ) +
  labs(
    x = "EZH2 cytoplsm",
    y = "CBX2 cytoplsm"
  ) + facet_wrap(~sample_type,scales = "free",labeller=facet_labeller) -> p1;p1
ggsave(filename = "EZH2-cytoplsm_CBX2-cytoplsm_immunehistochemistry_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 6,height = 3)
ggsave(filename = "EZH2-cytoplsm_CBX2-cytoplsm_immunehistochemistry_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 6,height = 3)


  
broom::tidy(
  cor.test(immune_histone.T$CBX2_cytoplsm,immune_histone.T$EZH2_karyon,method = "kendall")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=1,y=1.8,sample_type="1Tumor",n=nrow(immune_histone.T), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2_karyon.T
broom::tidy(
  cor.test(immune_histone.N$CBX2_cytoplsm,immune_histone.N$EZH2_karyon,method = "kendall")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=0.41,y=0.8,sample_type="2Normal",n=nrow(immune_histone.N), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2_karyon.N

rbind(CBX2_EZH2_karyon.T,CBX2_EZH2_karyon.N) %>%
  dplyr::arrange(sample_type) ->CBX2_EZH2_karyon

immune_histone %>%
  ggplot(aes(x=EZH2_karyon,y=CBX2_cytoplsm)) +
  geom_point(aes(color = sample_type)) +
  geom_smooth(se = FALSE, fullrange=TRUE, color = "#039BE5") +
  facet_wrap(~sample_type,scales = "free",labeller = facet_labeller) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  labs(
    x = "EZH2 karyon",
    y = "CBX2 cytoplsm"
  ) +
  scale_color_manual(
    values = c("#EE6363","#00C5CD"),
    labels = CBX2_EZH2_karyon$label
  ) -> p2;p2
ggsave(filename = "EZH2-karyo_CBX2-cytoplsm_immunehistochemistry_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 6,height = 3)
ggsave(filename = "EZH2-karyo_CBX2-cytoplsm_immunehistochemistry_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 6,height = 3)


library(grid)
library(gridExtra)
library(scales)
p2 + theme(axis.title.y = element_blank()) ->p2
# E Zhou -----
pdf(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.pdf"),width = 8,height = 3)
tiff(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.tiff"),width = 8,height = 3,units ="in",res=400)

# Home ------
pdf(file.path("D:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.pdf"),width = 8,height = 3)
tiff(file.path("D:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2_CBX2_immunehistochemistry_correlation.tiff"),width = 8,height = 3,units ="in",res=400)

grid.arrange(p1, p2, ncol=2,nrow=1,widths=c(1,1), heights=c(1))
dev.off()

# For CBX2 and EZH2 karyon
immune_histone.T %>%
  tidyr::drop_na() -> immune_histone.T.noNA
broom::tidy(
  cor.test(immune_histone.T.noNA$CBX2_karyon,immune_histone.T.noNA$EZH2_karyon,method = "kendall")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr"))%>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=1,y=1.8,sample_type="1Tumor",n=nrow(immune_histone.T.noNA), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2_karyon.T
immune_histone.N  %>%
  dplyr::select(-stage)%>%
  tidyr::drop_na()-> immune_histone.N.noNA
broom::tidy(
  cor.test(immune_histone.N.noNA$CBX2_karyon,immune_histone.N.noNA$EZH2_karyon,method = "kendall")) %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=0.41,y=0.8,sample_type="2Normal",n=nrow(immune_histone.N.noNA), estimate = signif(estimate,2)) %>%
  dplyr::mutate(label=purrr::map2(
    .x=p.value,
    .y=estimate,
    .z=n,
    .f=function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  )) ->CBX2_EZH2_karyon.N

rbind(CBX2_EZH2_karyon.T,CBX2_EZH2_karyon.N) ->CBX2_EZH2_karyon

immune_histone %>%
  dplyr::select(-stage) %>%
  tidyr::drop_na() %>%
  ggplot(aes(x=EZH2_karyon,y=CBX2_karyon)) +
  geom_point(aes(color = sample_type)) +
  geom_smooth(se = FALSE, fullrange=TRUE, color = "#039BE5") +
  facet_wrap(~sample_type,scales = "free",labeller = facet_labeller) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  labs(
    x = "EZH2 karyon",
    y = "CBX2 karyon"
  ) +
  scale_color_manual(
    values = c("#EE6363","#00C5CD"),
    labels = CBX2_EZH2_karyon$label
  ) -> p2;p2
ggsave(filename = "EZH2-karyo_CBX2-karyo_immunehistochemistry_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 6,height = 3)
ggsave(filename = "EZH2-karyo_CBX2-karyo_immunehistochemistry_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 6,height = 3)


#########################heatmap
rownames(immune_histone) <-immune_histone$sample_type %>% paste(c(1:nrow(immune_histone)),sep=".")

immune_histone %>%
  dplyr::select(CBX2_cytoplsm,EZH2_karyon,EZH2_cytoplsm) %>%
  pheatmap::pheatmap(cluster_cols = F)

##########################stage
oneway.test(EZH2_karyon~stage,data = immune_histone)
boxplot(EZH2_karyon~stage,data = immune_histone)

library(ggplot2)
comp_list <- list(c("Normal", "stage II"), c("stage II", "stage III"))
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
  dplyr::select(-c(EZH2_cytoplsm,EZH2_karyon,CBX2_cytoplsm,CBX2_karyon)) %>%
  tidyr::gather(-c("sample","sample_type","stage"),key="group",value="PositiveRateXStainingIntensity") %>%
  dplyr::mutate(sample_type=as.character(sample_type)) %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="Normal","1","2")) %>%
  dplyr::mutate(group=sub(pattern = "_",replacement = " ",group)) %>%
  dplyr::arrange(sample_type) %>%
  dplyr::group_by(group)->df
b <- runif(nrow(df), -0.2, 0.2)
comp_list <- list(c("2","1"))
df %>%
  ggpubr::ggboxplot(x = "sample_type", y = "PositiveRateXStainingIntensity",
                    color = "sample_type", pallete = "npg",add = "jitter",
                    facet.by = "group") +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("1","2"),
                   labels=c("Normal","Tumor")) +
  scale_color_manual(
    values = c("#00C5CD","#EE6363")
  )+
  theme(
    axis.title.x = element_blank()
  ) +
  ylab("PositiveRate X StainingIntensity") +
  scale_y_continuous(limits = c(0,2.5)) +
  # ggpubr::stat_compare_means(label.y = 2.3,method = "wilcox.test",label = "p.format") +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(2.5),label = "p.signif") ->p;p
ggsave(filename = "immune_histone.tiff", plot = p, device = "tiff", path = "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",
       width = 4, height = 3)
ggsave(filename = "immune_histone.pdf", plot = p, device = "pdf", path = "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",
       width = 4, height = 3)

save.image(file = file.path(result_path,"immune_histone.rdata"))
