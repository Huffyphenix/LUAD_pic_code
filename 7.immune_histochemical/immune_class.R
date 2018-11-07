.libPaths("E:/library")
library(magrittr)
library(ggplot2)
# HOME -----
data_path<-"Z:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"

# E Zhou -----
data_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"

# HUST ----
data_path<-"G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/data"
result_path<-"G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/芯片-免疫组化/result"


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
  dplyr::mutate(class = ifelse(all_cell_score>=6, "High", class)) -> immune_class

## for tumor samples
immune_class %>%
  dplyr::filter(sample_type == "T") %>%
  dplyr::select(sample,Gene,all_cell_score) %>%
  dplyr::arrange(Gene) %>%
  tidyr::spread(key = "Gene", value = c("all_cell_score")) -> immune_score.T

immune_class %>%
  dplyr::filter(sample_type == "T") %>%
  dplyr::arrange(Gene) %>%
  dplyr::mutate(class_n = ifelse(class == "High", "1_High", "2_middle")) %>%
  dplyr::mutate(class_n = ifelse(class == "Low", "3_Low", class_n)) %>%
  dplyr::select(sample,Gene,class_n) %>%
  tidyr::spread(key = "Gene", value = c("class_n")) -> immune_class.T

hospital_names <- list(
  "1_High" = "High",
  "2_middle" = "Middle", 
  "3_Low" = "Low"
)
hospital_labeller <- function(variable,value){
  return(hospital_names[value])
}

immune_score.T %>%
  dplyr::inner_join(immune_class.T, by = "sample") -> immune_class_score.T
immune_class_score.T %>%
  ggplot(aes(x=CBX2.x, y=EZH2.x)) +
  geom_jitter() +
  facet_grid(CBX2.y ~ EZH2.y, scales = "free")

immune_class_score.T %>%
  dplyr::filter(CBX2.y == "2_middle") %>%
  dplyr::filter(EZH2.y == "2_middle") -> immune_class_score.T.all_middle

#### 卡方检验
immune_class_score.T %>%
  dplyr::group_by(CBX2.y, EZH2.y) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(CBX2.y,EZH2.y,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=EZH2.y,value= n) -> immune_class_score.T.conti_table

xtabs(~CBX2.y+EZH2.y, data=immune_class_score.T) -> mytable
prop.table(mytable)
chisq.test(mytable,correct = T)


immune_class_score.T %>%
  dplyr::mutate(CBX2.y = ifelse(CBX2.y=="2_middle", CBX2.y, "H/L")) %>%
  dplyr::mutate(EZH2.y = ifelse(EZH2.y=="2_middle", CBX2.y, "H/L")) -> immune_class_score.T.2group

xtabs(~CBX2.y+EZH2.y, data=immune_class_score.T.2group) -> mytable
prop.table(mytable)
chisq.test(mytable,correct = T)

#### 相关性分析，增加随机扰动
set.seed(5000) # 设定种子
y <-rnorm(43,sd=0.05) 

immune_class.T %>%
  dplyr::inner_join(immune_class_score.T.all_middle,by="sample") %>%
  dplyr::mutate(CBX2.x=CBX2.x+y,EZH2.x=EZH2.x+y) -> immune_histone_raw.T.all_middle
broom::tidy(
  cor.test(immune_histone_raw.T.all_middle$CBX2.x,immune_histone_raw.T.all_middle$EZH2.x,method = "spearman")
)

## for normal samples 
immune_class %>%
  dplyr::filter(sample_type == "TA") %>%
  dplyr::select(sample,Gene,all_cell_score) %>%
  dplyr::arrange(Gene) %>%
  tidyr::spread(key = "Gene", value = c("all_cell_score")) -> immune_score.TA

immune_class %>%
  dplyr::filter(sample_type == "TA") %>%
  dplyr::arrange(Gene) %>%
  dplyr::mutate(class_n = ifelse(class == "High", "1_High", "2_middle")) %>%
  dplyr::mutate(class_n = ifelse(class == "Low", "3_Low", class_n)) %>%
  dplyr::select(sample,Gene,class_n) %>%
  tidyr::spread(key = "Gene", value = c("class_n")) -> immune_class.TA

immune_score.TA %>%
  dplyr::inner_join(immune_class.TA, by = "sample") -> immune_class_score.TA
immune_class_score.TA %>%
  ggplot(aes(x=CBX2.x, y=EZH2.x)) +
  geom_point() +
  geom_jitter(width=0.5,height=0.5) +
  facet_grid(CBX2.y ~ EZH2.y, scales = "free")


xtabs(~CBX2.y+EZH2.y, data=immune_class_score.TA) -> mytable
prop.table(mytable)
chisq.test(mytable,correct = T)



# do correlation for all middle samples -----------------------------------
# load data
# This data are raw data with staining * positive rate
immune_histone_raw  <- read.table(file.path(data_path,"immune_histone.txt"),sep = "\t",header = T) %>%
  dplyr::mutate(stage = as.character(stage)) %>%
  dplyr::mutate(stage = ifelse(is.na(stage),"N",stage)) %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="T","Tumor","Normal")) %>%
  dplyr::mutate(EZH2_cytoplsm=ifelse(is.na(EZH2_cytoplsm),0,EZH2_cytoplsm)) %>%
  dplyr::mutate(EZH2_karyon=ifelse(is.na(EZH2_karyon),0,EZH2_karyon)) %>%
  dplyr::mutate(CBX2_cytoplsm=ifelse(is.na(CBX2_cytoplsm),0,CBX2_cytoplsm)) %>%
  dplyr::mutate(CBX2_karyon=ifelse(is.na(CBX2_karyon),0,CBX2_karyon)) %>%
  dplyr::mutate(CBX2_mean = (CBX2_cytoplsm+CBX2_karyon)/2) %>%
  dplyr::mutate(EZH2_mean = (EZH2_karyon+EZH2_cytoplsm)/2) %>%
  dplyr::mutate(CBX2_max = ifelse(CBX2_cytoplsm > CBX2_karyon, CBX2_cytoplsm, CBX2_karyon)) %>%
  dplyr::mutate(EZH2_max = ifelse(EZH2_cytoplsm > EZH2_karyon, EZH2_cytoplsm, EZH2_karyon))



# Preleminary test to check the test assumptions
shapiro.test(immune_histone_raw.T.all_middle$CBX2_mean) # p-value < 0.05, don't follow a normal distribution.

broom::tidy(
  cor.test(immune_histone_raw.T.all_middle$CBX2_mean,immune_histone_raw.T.all_middle$EZH2_mean,method = "pearson")
)

ggplot(immune_histone_raw.T.all_middle,aes(x=CBX2_max,y=EZH2_max)) +
  geom_point() +
  geom_jitter()
