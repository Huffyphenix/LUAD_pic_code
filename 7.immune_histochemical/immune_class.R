.libPaths("E:/library")
library(magrittr)
library(ggplot2)
# HOME -----
data_path<-"Z:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"

# E Zhou -----
data_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2"

# HUST ----
data_path<-"G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2"


# load data ---------------------------------------------------------------
immune_histone  <- read.table(file.path(data_path,"immune_histone_data_for_classify.txt"),sep = "\t",header = T)

immune_histone %>% # negative and NA all represent by 0.
  dplyr::mutate(positiveRate_c_score = ifelse(positiveRate_c<0.25,1,2)) %>%
  dplyr::mutate(positiveRate_c_score = ifelse(positiveRate_c>=0.5,3,positiveRate_c_score)) %>%
  dplyr::mutate(positiveRate_k_score = ifelse(positiveRate_k<0.25,1,2)) %>%
  dplyr::mutate(positiveRate_k_score = ifelse(positiveRate_k>=0.5,3,positiveRate_k_score)) %>%
  dplyr::mutate(k = positiveRate_k_score * staningIntensity_k) %>%
  dplyr::mutate(c = positiveRate_c_score * staningIntensity_c) %>%
  dplyr::mutate(all_cell_score = ifelse(k>c,k,c)) %>%
  dplyr::mutate(class = ifelse(all_cell_score<3, "Low", "Middle")) %>%
  dplyr::mutate(class = ifelse(all_cell_score>=6, "High", class)) %>%
  dplyr::mutate(sample_type = ifelse(sample_type == "T", "Tumor", "Normal")) -> immune_class

# Score for positive rate
# From reference: https://www.spandidos-publications.com/10.3892/mmr.2016.5714
# 
immune_histone %>% # negative and NA all represent by 0.
  dplyr::mutate(positiveRate_c_score = ifelse(positiveRate_c<0.05,0,1)) %>%
  dplyr::mutate(positiveRate_c_score = ifelse(positiveRate_c>=0.25 & positiveRate_c <0.5,2,positiveRate_c_score)) %>%
  dplyr::mutate(positiveRate_c_score = ifelse(positiveRate_c>=0.5,3,positiveRate_c_score)) %>%
  dplyr::mutate(positiveRate_k_score = ifelse(positiveRate_k<0.05,0,1)) %>%
  dplyr::mutate(positiveRate_k_score = ifelse(positiveRate_k>=0.25 & positiveRate_k <0.5,2,positiveRate_k)) %>%
  dplyr::mutate(positiveRate_k_score = ifelse(positiveRate_k>=0.5, 3,positiveRate_k)) %>%
  dplyr::mutate(k = positiveRate_k_score * staningIntensity_k) %>%
  dplyr::mutate(c = positiveRate_c_score * staningIntensity_c) %>%
  dplyr::mutate(all_cell_score = ifelse(k>c,k,c)) %>%
  dplyr::mutate(class = ifelse(all_cell_score<3, "Low", "Middle")) %>%
  dplyr::mutate(class = ifelse(all_cell_score>=6, "High", class)) %>%
  dplyr::mutate(sample_type = ifelse(sample_type == "T", "Tumor", "Normal")) -> immune_class
######################################################################
## annotate the code seprerate calculate the tumor and normal samples
## do class, correaltion and generate plot for all samples togather
#####################################################################
## for tumor samples
# immune_class %>%
#   dplyr::filter(sample_type == "T") %>%
#   dplyr::select(sample, Gene,all_cell_score) %>%
#   dplyr::arrange(Gene) %>%
#   tidyr::spread(key = "Gene", value = c("all_cell_score")) -> immune_score.T
# 
# immune_class %>%
#   dplyr::filter(sample_type == "T") %>%
#   dplyr::arrange(Gene) %>%
#   dplyr::mutate(class_n = ifelse(class == "High", "1_High", "2_middle")) %>%
#   dplyr::mutate(class_n = ifelse(class == "Low", "3_Low", class_n)) %>%
#   dplyr::select(sample,Gene,class_n) %>%
#   tidyr::spread(key = "Gene", value = c("class_n")) -> immune_class.T

# ggplot
# hospital_names <- list(
#   "1_High" = "High",
#   "2_middle" = "Middle", 
#   "3_Low" = "Low"
# )
# hospital_labeller <- function(variable,value){
#   return(hospital_names[value])
# }
# 
# immune_score.T %>%
#   dplyr::inner_join(immune_class.T, by = "sample") -> immune_class_score.T
# 
# ### 1 
# immune_class_score.T %>%
#   dplyr::mutate(group_combine = ifelse(CBX2.y=="2_middle" & EZH2.y=="2_middle", "1", "2")) %>%
#   dplyr::group_by(CBX2.y,EZH2.y) %>%
#   dplyr::mutate(n=paste("n = ",n(), sep = "")) %>%
#   dplyr::ungroup() %>%
#   ggplot(aes(x=CBX2.x, y=EZH2.x)) +
#   geom_jitter(aes(color = group_combine),width = 1, height = 1) +
#   facet_grid(CBX2.y ~ EZH2.y, scales = "free",labeller = hospital_labeller) +
#   geom_text(aes(label = n)) +
#   scale_color_manual(values = c("#FF3030", "#050505")) +
#   theme(
#     plot.background = element_rect(fill = "white", colour = "black"),
#     panel.background = element_rect(fill = "white", colour = "black"),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black", size = 0.5),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 12),
#     axis.title = element_text(size = 12, color = "black"),
#     axis.text = element_text(size = 12, colour = "black"),
#     legend.position = "none"
#   ) +
#   ylab("EZH2 protein level") +
#   xlab("CBX2 protein level")
# ggsave(file.path(result_path,"immune_histone_T_distribution.pdf"),device = "pdf",height = 3,width = 4)
# ggsave(file.path(result_path,"immune_histone_T_distribution.tiff"),device = "tiff",height = 3,width = 4)
# 
# immune_class_score.T %>%
#   dplyr::filter(CBX2.y == "2_middle") %>%
#   dplyr::filter(EZH2.y == "2_middle") -> immune_class_score.T.all_middle
# 
# #### 卡方检验，两个group
# immune_class_score.T %>%
#   dplyr::group_by(CBX2.y, EZH2.y) %>%
#   dplyr::mutate(n = n()) %>%
#   dplyr::select(CBX2.y,EZH2.y,n) %>%
#   unique() %>%
#   dplyr::ungroup() %>%
#   tidyr::spread(key=EZH2.y,value= n) -> immune_class_score.T.conti_table
# 
# xtabs(~CBX2.y+EZH2.y, data=immune_class_score.T) -> mytable
# prop.table(mytable)
# chisq.test(mytable,correct = T)
# 
# 
# immune_class_score.T %>%
#   dplyr::mutate(CBX2.y = ifelse(CBX2.y=="2_middle", CBX2.y, "H/L")) %>%
#   dplyr::mutate(EZH2.y = ifelse(EZH2.y=="2_middle", CBX2.y, "H/L")) -> immune_class_score.T.2group
# 
# xtabs(~CBX2.y+EZH2.y, data=immune_class_score.T.2group) -> mytable
# prop.table(mytable)
# chisq.test(mytable,correct = T)
# 
# #### 相关性分析，增加随机扰动
# 
# # For all middle samples
# set.seed(5000) # 设定种子
# y <-rnorm(43,sd=0.01) 
# 
# immune_class.T %>%
#   dplyr::inner_join(immune_class_score.T.all_middle,by="sample") %>%
#   dplyr::mutate(CBX2.x=CBX2.x+y,EZH2.x=EZH2.x+y) -> immune_histone_raw.T.all_middle
# broom::tidy(
#   cor.test(immune_histone_raw.T.all_middle$CBX2.x,immune_histone_raw.T.all_middle$EZH2.x,method = "spearman")
# )
# 
# immune_histone_raw.T.all_middle %>%
#   ggplot(aes(x=CBX2.x, y=EZH2.x)) +
#   geom_jitter() +
#   scale_color_manual(values = c("#FF3030", "#050505")) +
#   theme(
#     plot.background = element_rect(fill = "white", colour = "black"),
#     panel.background = element_rect(fill = "white", colour = "black"),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black", size = 0.5),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 12),
#     axis.title = element_text(size = 12, color = "black"),
#     axis.text = element_text(size = 12, colour = "black"),
#     legend.position = "none"
#   ) +
#   ylab("EZH2 protein level") +
#   xlab("CBX2 protein level")
# # For all samples
# set.seed(5000) # 设定种子
# y <-rnorm(75,sd=0.01) 
# immune_class_score.T %>%
#   dplyr::mutate(CBX2.x=CBX2.x+y,EZH2.x=EZH2.x+y) ->immune_class_score.T.corrlate
# broom::tidy(
#   cor.test(immune_class_score.T.corrlate$CBX2.x,immune_class_score.T.corrlate$EZH2.x,method = "spearman")
# )
# immune_class_score.T %>%
#   dplyr::mutate(group_combine = ifelse(CBX2.y=="2_middle" & EZH2.y=="2_middle", "1", "2")) %>%
#   ggplot(aes(x=CBX2.x, y=EZH2.x)) +
#   geom_jitter(height = 0.1, width = 0.1,aes(color = group_combine)) +
#   scale_color_manual(values = c("#FF3030", "#050505")) +
#   geom_smooth(method = "lm") +
#   theme(
#     panel.background = element_rect(fill = "white"),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black", size = 0.5),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 12),
#     axis.title = element_text(size = 12, color = "black"),
#     axis.text = element_text(size = 12, colour = "black"),
#     legend.position = "none"
#   ) +
#   ylab("EZH2 protein level") +
#   xlab("CBX2 protein level")
# 
# 
# 
# # for normal samples  -----------------------------------------------------
# ## for normal samples 
# immune_class %>%
#   dplyr::filter(sample_type == "TA") %>%
#   dplyr::select(sample,Gene,all_cell_score) %>%
#   dplyr::arrange(Gene) %>%
#   tidyr::spread(key = "Gene", value = c("all_cell_score")) -> immune_score.TA
# 
# 
# 
# 
# immune_class %>%
#   dplyr::filter(sample_type == "TA") %>%
#   dplyr::arrange(Gene) %>%
#   dplyr::mutate(class_n = ifelse(class == "High", "1_High", "2_middle")) %>%
#   dplyr::mutate(class_n = ifelse(class == "Low", "3_Low", class_n)) %>%
#   dplyr::select(sample,Gene,class_n) %>%
#   tidyr::spread(key = "Gene", value = c("class_n")) -> immune_class.TA
# 
# immune_score.TA %>%
#   dplyr::inner_join(immune_class.TA, by = "sample") -> immune_class_score.TA
# immune_class_score.TA %>%
#   ggplot(aes(x=CBX2.x, y=EZH2.x)) +
#   geom_point() +
#   geom_jitter(width=0.5,height=0.5) +
#   facet_grid(CBX2.y ~ EZH2.y, scales = "free")
# 
# 
# xtabs(~CBX2.y+EZH2.y, data=immune_class_score.TA) -> mytable
# prop.table(mytable)
# chisq.test(mytable,correct = T)



# For tumor and normal -----------------------------------
immune_class %>%
  dplyr::select(sample,sample_type, Gene,all_cell_score) %>%
  dplyr::arrange(Gene) %>%
  tidyr::spread(key = "Gene", value = c("all_cell_score")) -> immune_score.all

immune_class %>%
  dplyr::arrange(Gene) %>%
  dplyr::mutate(class_n = ifelse(class == "High", "1_High", "2_middle")) %>%
  dplyr::mutate(class_n = ifelse(class == "Low", "3_Low", class_n)) %>%
  dplyr::select(sample,Gene,class_n) %>%
  tidyr::spread(key = "Gene", value = c("class_n")) -> immune_class.all

# ggplot
hospital_names <- list(
  "1_High" = "High",
  "2_middle" = "Middle", 
  "3_Low" = "Low",
  "Normal" = "Normal",
  "Tumor" = "Tumor"
)
hospital_labeller <- function(variable,value){
  return(hospital_names[value])
}

immune_score.all %>%
  dplyr::inner_join(immune_class.all, by = "sample") -> immune_class_score.all

### plot to see the distribution of CBX2 and EZH2 expression
immune_class_score.all %>%
  dplyr::mutate(group_combine = ifelse(CBX2.y=="2_middle" & EZH2.y=="2_middle" & sample_type == "Tumor", "1", "2")) %>%
  dplyr::mutate(group_combine = ifelse(CBX2.y=="3_Low" & EZH2.y=="3_Low" & sample_type == "Normal", "1",group_combine)) %>%
  dplyr::group_by(sample_type,CBX2.y,EZH2.y) %>%
  dplyr::mutate(n=paste("n = ",n(), sep = "")) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x=CBX2.x, y=EZH2.x)) +
  geom_jitter(aes(color = group_combine, shape = sample_type),width = 1, height = 1, size = 0.5) +
  # facet_wrap(~sample_type) +
  facet_grid( EZH2.y ~ sample_type+CBX2.y, labeller = hospital_labeller ) +
  geom_text(aes(label = n)) +
  scale_color_manual(values = c("#FF3030", "#050505")) +
  scale_shape_manual(values = c(1,2)) +
  theme(
    plot.background = element_rect(fill = "white", colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = "none"
  ) +
  # xlim(0,6) +
  ylab("EZH2 protein level") +
  xlab("CBX2 protein level")
ggsave(file.path(result_path,"immune_histone_T-N_distribution-1.pdf"),device = "pdf",height = 3,width = 5)
ggsave(file.path(result_path,"immune_histone_T-N_distribution.tiff"),device = "tiff",height = 4,width = 5)


human_read <- function(.x){
  if (.x > 0.1) {
    .x %>% signif(digits = 2) %>% toString()
  } else if (.x < 0.1 && .x > 0.001 ) {
    .x %>% signif(digits = 1) %>% toString()
  } else {
    .x %>% format(digits = 2, scientific = TRUE)
  }
}

# do correlation For all samples
set.seed(5000) # 设定种子, make the results repeatable
y <-rnorm(150,sd=0.01) # add a Random disturbance for all samples
immune_class_score.all %>%
  dplyr::mutate(CBX2.x=CBX2.x+y,EZH2.x=EZH2.x+y) ->immune_class_score.all.corrlate
immune_class_score.all.corrlate %>%
  dplyr::group_by(sample_type) %>%
  dplyr::do(
    cor = broom::tidy(
      cor.test(.$CBX2.x,.$EZH2.x,method = "spearman")
    ), 
    n = nrow(.)
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value,human_read)) %>%
  dplyr::mutate(x=1,y=1.8, estimate = signif(estimate,2)) %>%
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
  )) -> cor_p_label

immune_class_score.all %>%
  dplyr::mutate(group_combine = ifelse(CBX2.y=="2_middle" & EZH2.y=="2_middle" & sample_type == "Tumor", "1", "2")) %>%
  dplyr::mutate(group_combine = ifelse(CBX2.y=="3_Low" & EZH2.y=="3_Low" & sample_type == "Normal", "1",group_combine)) %>%
  ggplot(aes(x=CBX2.x, y=EZH2.x)) +
  geom_jitter(height = 0.1, width = 0.1,aes(color = sample_type, shape = sample_type)) +
  scale_shape_manual(values = c(1,2),
                     labels = cor_p_label$label) +
  facet_wrap(~ sample_type) +
  scale_color_manual(values = c("#00BFFF","#EE6363")) +
  geom_smooth(method = "lm") +
  # geom_vline(xintercept = c(2.8,5.9)) +
  # geom_hline(yintercept = c(5.9,2.8))+
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = c(0.2,0.8),
    legend.background = element_blank()
  ) +
  ylab("EZH2 protein level") +
  xlab("CBX2 protein level")

ggsave(file.path(result_path,"immune_histone_T-N_correlation.pdf"),device = "pdf",height = 3,width = 3)
ggsave(file.path(result_path,"immune_histone_T-N_correlation.tiff"),device = "tiff",height = 4,width = 5)

# DE analysis between tumor and normal samples
comp_list <- list(c("2Tumor","1Normal"))
immune_class %>%
  dplyr::select(sample,sample_type,Gene,all_cell_score,class) %>%
  dplyr::mutate(sample_type = ifelse(sample_type=="Tumor","2Tumor","1Normal")) %>%
  dplyr::arrange(sample_type) %>%
  ggpubr::ggboxplot(x = "sample_type", y = "all_cell_score",
                    color = "sample_type", palette = "npg" #add = "jitter",
  ) +
  facet_wrap(~Gene, strip.position = "bottom") +
  geom_jitter(aes(color=sample_type)) +
  scale_color_manual(
    values = c("#1E90FF","#EE6363")
  )+
  ylim(0,7)+
  ylab("Protein level") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white",colour = "white"),
        strip.text = element_text(size = 12)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(6.5),label = "p.signif") +
  scale_x_discrete(breaks = c("1Normal","2Tumor"),
                   labels = c("Normal","Tumor"),
                   expand = c(0.2,0.2))
ggsave(file.path(result_path,"immune_histone_T-N_DE.pdf"),device = "pdf",height = 3,width = 5)
ggsave(file.path(result_path,"immune_histone_T-N_DE.tiff"),device = "tiff",height = 3,width = 5)
