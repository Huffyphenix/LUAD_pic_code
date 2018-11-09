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

immune_histone  <- read.table(file.path(data_path,"immune_histone.txt"),sep = "\t",header = T) %>%
  dplyr::mutate(stage = as.character(stage)) %>%
  dplyr::mutate(stage = ifelse(is.na(stage),"N",stage)) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(sample_type=ifelse(sample_type=="T","Tumor","Normal")) %>%
  # dplyr::mutate(EZH2_cytoplsm=ifelse(is.na(EZH2_cytoplsm),0,EZH2_cytoplsm)) %>%
  # dplyr::mutate(EZH2_karyon=ifelse(is.na(EZH2_karyon),0,EZH2_karyon)) %>%
  # dplyr::mutate(CBX2_cytoplsm=ifelse(is.na(CBX2_cytoplsm),0,CBX2_cytoplsm)) %>%
  # dplyr::mutate(CBX2_karyon=ifelse(is.na(CBX2_karyon),0,CBX2_karyon)) %>%
  dplyr::mutate(CBX2_mean = (CBX2_cytoplsm+CBX2_karyon)/2) %>%
  dplyr::mutate(EZH2_mean = (EZH2_karyon+EZH2_cytoplsm)/2) %>%
  dplyr::mutate(CBX2_max = ifelse(CBX2_cytoplsm > CBX2_karyon, CBX2_cytoplsm, CBX2_karyon)) %>%
  dplyr::mutate(EZH2_max = ifelse(EZH2_cytoplsm > EZH2_karyon, EZH2_cytoplsm, EZH2_karyon))
  
# Preleminary test to check the test assumptions
shapiro.test(immune_histone$EZH2_mean) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(immune_histone$CBX2_mean) # p-value < 0.05, don't follow a normal distribution.


# do correlation for tumor and normal samples, respectively
## data precessing
immune_histone %>%
  dplyr::filter(sample_type=="Tumor") -> immune_histone.T
immune_histone %>%
  dplyr::filter(sample_type=="Normal") -> immune_histone.N

## function to get scientific numeric 
human_read <- function(.x){
  if (.x > 0.1) {
    .x %>% signif(digits = 2) %>% toString()
  } else if (.x < 0.1 && .x > 0.001 ) {
    .x %>% signif(digits = 1) %>% toString()
  } else {
    .x %>% format(digits = 2, scientific = TRUE)
  }
}

## do

broom::tidy(
  cor.test(immune_histone.T$CBX2_mean,immune_histone.T$EZH2_mean,method = "pearson")) %>%
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
  cor.test(immune_histone.N$CBX2_max,immune_histone.N$EZH2_max,method = "kendall")) %>%
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
  )) ->CBX2_EZH2.N
rbind(CBX2_EZH2.T,CBX2_EZH2.N) %>%
  dplyr::as.tbl() ->CBX2_EZH2;CBX2_EZH2


facet_names <- list(
  '1Tumor'="Tumor",
  '2Normal'="Normal"
)
facet_labeller <- function(variable,value){
  return(facet_names[value])
}


immune_histone %>%
  ggplot(aes(x=EZH2_mean,y=CBX2_mean)) +
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
    labels = CBX2_EZH2$label
  ) +
  labs(
    x = "EZH2 cytoplsm",
    y = "CBX2 cytoplsm"
  ) + facet_wrap(~sample_type,scales = "free",labeller=facet_labeller) -> p1;p1
