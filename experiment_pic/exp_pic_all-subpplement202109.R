library(magrittr)
library(plyr)
library(ggplot2)
data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/03.LUAD/02.实验图片/整理实验图/2021-09补实验/process_data"
result_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/03.LUAD/02.实验图片/整理实验图/2021-09补实验/pic_by_R"
### function -------------------------------------
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
##################################################
### CBX2 and EZH2 siRNA cell variability in A549

viability <- readr::read_tsv(file.path(data_path,"MTT.txt")) %>%
  tidyr::gather(-h,-group,key = "rep",value="viability")

data_summary(viability,varname = "viability",groupnames = c("group","h")) -> viability_summary

# get p signif
viability %>%
  tidyr::spread(key="group",value="viability") %>%
  dplyr::group_by(h) %>%
  tidyr::nest() %>%
  dplyr::filter(h!="0") %>%
  dplyr::mutate(pval=purrr::map(data,.f=function(.x){
    
    .p <- data.frame(group1=c("siControl+DMSO","siCBX2+DMSO","siCBX2+TAZ","siControl+TAZ"),
                        group2=c("siCBX2+DMSO","siCBX2+TAZ","siControl+TAZ","siControl+DMSO"),p=NA)
    for(i in 1:nrow(.p)){
      t.test(.x[,.p[i,"group1"]],.x[,.p[i,"group2"]]) %>% 
        broom::tidy() %>% .[1,5] -> .p[i,"p"]
    }
    .p
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> viability_ttest

viability_ttest %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
    dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe))-> viability_ttest_plabel
viability_ttest_plabel %>%
  readr::write_tsv(file.path(result_path,"viability_ttest_plabel.txt"))
# ggplot
viability_summary<-within(viability_summary,group <- factor(group,levels = c("siControl+DMSO","siCBX2+DMSO","siControl+TAZ","siCBX2+TAZ")))
with(viability_summary,levels(group))
ggplot(viability_summary, aes(x=group, y=viability, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~h,strip.position = "bottom",nrow = 1) +
  # geom_text(aes(y = viability+sd+1, label=p_labe),size=5) +
  theme_classic() +
  scale_fill_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend.key.width=unit(0.15,"inches"),  # legend size
    # legend.key.height=unit(0.15,"inches"),
    strip.text = element_text(size = 15),
    strip.background = element_rect(colour = "white")
  ) +
  ylab("OD 490 nm") +
  xlab("Time (hours)")+
  labs(title = "A549 cell viability")

ggsave(file.path(result_path,"viability_siCBX-TAZ_A549_bar.pdf"),width = 7,height = 3,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-TAZ_A549_bar.tiff"),width = 7,height = 3,device = "tiff")

## broken line
viability_summary %>%
  dplyr::mutate(`Time (hours)`=as.character(h)) %>%
  ggplot(aes(x=`Time (hours)`, y=viability, color=group)) +
  geom_line(aes(group=group)) +
  geom_point() +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.1) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_color_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(color = "black", size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width=unit(0.15,"inches"),  # legend size
    legend.key.height=unit(0.15,"inches")
  )+
  ylab("OD 490 nm") +
  xlab("Time (hours)")

ggsave(file.path(result_path,"viability_siCBX-TAZ_A549.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-TAZ_A549.tiff"),width = 4,height = 2,device = "pdf")


##################################################
### tranwell_supplement
tranwell_supplement <- readr::read_tsv(file.path(data_path,"transwell.txt"))  %>%
  tidyr::gather(-group,key="rep",value="value")

data_summary(tranwell_supplement,varname = "value",groupnames = c("group")) -> tranwell_summary

# get p signif
tranwell_supplement %>%
  tidyr::spread(key="group",value="value") %>%
  tidyr::nest() %>%
  dplyr::mutate(pval=purrr::map(data,.f=function(.x){
    
    .p <- data.frame(group1=c("siControl+DMSO","siCBX2+DMSO","siCBX2+TAZ","siControl+TAZ"),
                     group2=c("siCBX2+DMSO","siCBX2+TAZ","siControl+TAZ","siControl+DMSO"),p=NA)
    for(i in 1:nrow(.p)){
      t.test(.x[,.p[i,"group1"]],.x[,.p[i,"group2"]]) %>% 
        broom::tidy() %>% .[1,5] -> .p[i,"p"]
    }
    .p
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tranwell_ttest

tranwell_ttest %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe))-> tranwell_ttest_plabel
tranwell_ttest_plabel %>%
  readr::write_tsv(file.path(result_path,"tranwell_ttest_plabel.txt"))
# ggplot
tranwell_summary<-within(tranwell_summary,group <- factor(group,levels = c("siControl+DMSO","siCBX2+DMSO","siControl+TAZ","siCBX2+TAZ")))
with(tranwell_summary,levels(group))
tranwell_summary %>%
  ggplot(aes(x=group, y=value, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#87CEFA", "#458B74")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) 

ggsave(file.path(result_path,"transwell_siCB2-TAZ.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"transwell_siCB2-TAZ.tiff"),width = 4,height = 2,device = "tiff")
