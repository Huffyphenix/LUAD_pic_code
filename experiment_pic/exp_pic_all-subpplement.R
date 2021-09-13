.libPaths("E:/library")
.libPaths("F:/library")
.libPaths("C:/Users/94998/Documents/library")
library(ggplot2)
library(magrittr)
library(plyr)
# data path ---------------------------------------------------------------
# HUST 
data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/process_data"
result_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/Pic_by_R"

# E Zhou 
data_path <- "E:/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/process_data"
result_path <- "E:/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/Pic_by_R"

# xiaomi laptop 
data_path <- "C:/Users/94998/Documents/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/process_data"
result_path <- "C:/Users/94998/Documents/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/Pic_by_R"

# wust 
data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/process_data"
result_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/Pic_by_R"

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
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
### CBX2 and EZH2 siRNA EDU in A549 cell ----
EDU_A549 <- readr::read_tsv(file.path(data_path,"1.EDU_A549.txt"),col_names = F) %>%
  tidyr::gather(-X1,key = "group",value="Relative_mRNA_level") %>%
  dplyr::select(-group) %>%
  dplyr::rename("group" = "X1")

data_summary(EDU_A549,varname = "Relative_mRNA_level",groupnames = "group") -> EDU_A549_summary

#  bar plot
# get p signif
EDU_A549 %>%
  dplyr::mutate(x=c(rep(1,4),rep(2,4),rep(3,4))) %>%
  tidyr::spread(key="group",value="Relative_mRNA_level") %>%
  as.matrix() -> EDU_A549_for_ttest
edu_p <- data.frame(group=EDU_A549_summary$group,p=NA)
for(i in 3:5){
  t.test(EDU_A549_for_ttest[,2],EDU_A549_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> edu_p[i-1,2]
}

edu_p_2 <- tibble::tibble()
for(i in EDU_A549_summary$group){
  for (j in EDU_A549_summary$group) {
    if(i!=j){
      t.test(EDU_A549_for_ttest[,i],EDU_A549_for_ttest[,j]) %>% 
        broom::tidy() %>% .[1,5] -> .p
      tibble::tibble(group1=i,group2=j,p=.p$p.value) -> tmp
      rbind(edu_p_2,tmp)->edu_p_2
    }else{
      edu_p_2->edu_p_2
    }
  }
}

edu_p_2 %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  readr::write_tsv(file.path(result_path,"EDU_siCBX-EZH2_A549-Pvalue.tsv"))

EDU_A549_summary %>%
  dplyr::inner_join(edu_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="A549 cell")-> EDU_A549_summary

# ggplot
EDU_A549_summary<-within(EDU_A549_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(EDU_A549_summary,levels(group))
ggplot(EDU_A549_summary, aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~x,strip.position = "bottom") +
  geom_text(aes(y = Relative_mRNA_level+sd+1, label=p_labe),size=5) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab(paste("EdU positive", "cells(%)",sep = "\n"))

ggsave(file.path(result_path,"EDU_siCBX-EZH2_A549.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"EDU_siCBX-EZH2_A549.tiff"),width = 4,height = 2,device = "tiff")

### CBX2 and EZH2 siRNA EDU in H1299 cell ----
EDU_H1299 <- readr::read_tsv(file.path(data_path,"1.EDU_H1299.txt"),col_names = F) %>%
  tidyr::gather(-X1,key = "group",value="Relative_mRNA_level") %>%
  dplyr::select(-group) %>%
  dplyr::rename("group" = "X1")

data_summary(EDU_H1299,varname = "Relative_mRNA_level",groupnames = "group") -> EDU_H1299_summary

#  bar plot
# get p signif
EDU_H1299 %>%
  dplyr::mutate(x=c(rep(1,4),rep(2,4),rep(3,4))) %>%
  tidyr::spread(key="group",value="Relative_mRNA_level") %>%
  as.matrix() -> EDU_H1299_for_ttest
edu_p <- data.frame(group=EDU_H1299_summary$group,p=NA)
for(i in 3:5){
  t.test(EDU_H1299_for_ttest[,2],EDU_H1299_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> edu_p[i-1,2]
}

edu_p_2 <- tibble::tibble()
for(i in EDU_H1299_summary$group){
  for (j in EDU_H1299_summary$group) {
    if(i!=j){
      t.test(EDU_H1299_for_ttest[,i],EDU_H1299_for_ttest[,j]) %>% 
        broom::tidy() %>% .[1,5] -> .p
      tibble::tibble(group1=i,group2=j,p=.p$p.value) -> tmp
      rbind(edu_p_2,tmp)->edu_p_2
    }else{
      edu_p_2->edu_p_2
    }
  }
}

edu_p_2 %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  readr::write_tsv(file.path(result_path,"EDU_siCBX-EZH2_H1299-Pvalue.tsv"))

EDU_H1299_summary %>%
  dplyr::inner_join(edu_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="H1299 cell")-> EDU_H1299_summary

# ggplot
EDU_H1299_summary<-within(EDU_H1299_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(EDU_H1299_summary,levels(group))
ggplot(EDU_H1299_summary, aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~x,strip.position = "bottom") +
  geom_text(aes(y = Relative_mRNA_level+sd+1, label=p_labe),size=5) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab(paste("EdU positive", "cells(%)",sep = "\n"))

ggsave(file.path(result_path,"EDU_siCBX-EZH2_H1299.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"EDU_siCBX-EZH2_H1299.tiff"),width = 4,height = 2,device = "tiff")


##### combine EDU of A549 and H1299 into one ----
EDU_H1299_summary %>%
  rbind(EDU_A549_summary) %>%
  ggplot(aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ x, strip.position = "left",ncol=1) +
  geom_text(aes(y = Relative_mRNA_level+sd+1.5, label=p_labe), size = 5, angle = 90) +
  # ylim(0,1.5)+
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank() ,
    axis.ticks.y = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.x = element_text(size = 12, colour = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size=12),
    # legend.key.height = unit(0.15,"inches"),
    # legend.key.width = unit(0.15,"inches"),
    legend.background = element_blank()
  )  +
  ylab(paste("EdU positive", "cells(%)")) +
  # labs(title = "Edu") +
  coord_flip()
ggsave(file.path(result_path,"EDU_siCBX-EZH2-H1299-A549.pdf"),width = 6,height = 2,device = "pdf")
ggsave(file.path(result_path,"EDU_siCBX-EZH2-H1299-A549.tiff"),width = 6,height = 2,device = "tiff")


ggsave(file.path(result_path,"EDU_siCBX-EZH2-H1299-A549-flip.tiff"),width = 4,height = 4,device = "tiff")
##################################################
### CBX2 and EZH2 siRNA invasion in A549 cell ----
A549_invasion <- readr::read_tsv(file.path(data_path,"4.Invasion_A549.txt")) %>%
  tidyr::gather(key = "group",value="invasion") %>%
  dplyr::mutate(invasion = invasion*100)

data_summary(A549_invasion,varname = "invasion",groupnames = "group") -> A549_invasion_summary

#  bar plot
# get p signif
A549_invasion %>%
  dplyr::mutate(x=c(rep(c(1,2,3),4))) %>%
  tidyr::spread(key="group",value="invasion") %>%
  as.matrix() -> A549_invasion_for_ttest
A549_invasion_p <- data.frame(group=A549_invasion_summary$group,p=NA)
for(i in 3:5){
  t.test(A549_invasion_for_ttest[,2],A549_invasion_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> A549_invasion_p[i-1,2]
}
A549_invasion_p_2 <- tibble::tibble()
for(i in A549_invasion_summary$group){
  for (j in A549_invasion_summary$group) {
    if(i!=j){
      t.test(A549_invasion_for_ttest[,i],A549_invasion_for_ttest[,j]) %>% 
        broom::tidy() %>% .[1,5] -> .p
      tibble::tibble(group1=i,group2=j,p=.p$p.value) -> tmp
      rbind(A549_invasion_p_2,tmp)->A549_invasion_p_2
    }else{
      A549_invasion_p_2->A549_invasion_p_2
    }
  }
}

A549_invasion_p_2 %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  readr::write_tsv(file.path(result_path,"invasion_siCBX-EZH2-A549-Pvalue.tsv"))

A549_invasion_summary %>%
  dplyr::inner_join(A549_invasion_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="A549 cell")-> A549_invasion_summary

# ggplot
A549_invasion_summary<-within(A549_invasion_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(A549_invasion_summary,levels(group))
ggplot(A549_invasion_summary, aes(x=group, y=invasion, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=invasion-sd, ymax=invasion+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ x, strip.position = "bottom") +
  geom_text(aes(y = invasion+sd+5, label=p_labe), size = 5) +
  # ylim(0,1.5)+
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size=12),
    # legend.key.height = unit(0.15,"inches"),
    # legend.key.width = unit(0.15,"inches"),
    legend.background = element_blank()
  )  +
  ylab(paste("Cells invasion", "(% of Control)",sep="\n"))

ggsave(file.path(result_path,"invasion_siCBX-EZH2-A549.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"invasion_siCBX-EZH2-A549.tiff"),width = 4,height = 2,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA invasion in H1299 cell
H1299_invasion <- readr::read_tsv(file.path(data_path,"4.Invasion_H1299.txt")) %>%
  tidyr::gather(key = "group",value="invasion") %>%
  dplyr::mutate(invasion = invasion*100)

data_summary(H1299_invasion,varname = "invasion",groupnames = "group") -> H1299_invasion_summary

.#  bar plot
# get p signif
H1299_invasion %>%
  dplyr::mutate(x=c(rep(c(1,2,3),4))) %>%
  tidyr::spread(key="group",value="invasion") %>%
  as.matrix() -> H1299_invasion_for_ttest
H1299_invasion_p <- data.frame(group=H1299_invasion_summary$group,p=NA)
for(i in 3:5){
  t.test(H1299_invasion_for_ttest[,2],H1299_invasion_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> H1299_invasion_p[i-1,2]
}
H1299_invasion_p_2 <- tibble::tibble()
for(i in H1299_invasion_summary$group){
  for (j in H1299_invasion_summary$group) {
    if(i!=j){
      t.test(H1299_invasion_for_ttest[,i],H1299_invasion_for_ttest[,j]) %>% 
        broom::tidy() %>% .[1,5] -> .p
      tibble::tibble(group1=i,group2=j,p=.p$p.value) -> tmp
      rbind(H1299_invasion_p_2,tmp)->H1299_invasion_p_2
    }else{
      H1299_invasion_p_2->H1299_invasion_p_2
    }
  }
}

H1299_invasion_p_2 %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  readr::write_tsv(file.path(result_path,"invasion_siCBX-EZH2-H1299-Pvalue.tsv"))
H1299_invasion_summary %>%
  dplyr::inner_join(H1299_invasion_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="H1299 cell")-> H1299_invasion_summary

# ggplot
H1299_invasion_summary<-within(H1299_invasion_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(H1299_invasion_summary,levels(group))
ggplot(H1299_invasion_summary, aes(x=group, y=invasion, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=invasion-sd, ymax=invasion+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ x, strip.position = "bottom") +
  geom_text(aes(y = invasion+sd+5, label=p_labe), size = 5) +
  # ylim(0,1.5)+
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size=12),
    # legend.key.height = unit(0.15,"inches"),
    # legend.key.width = unit(0.15,"inches"),
    legend.background = element_blank()
  )  +
  ylab(paste("Cells invasion", "(% of Control)",sep="\n"))

ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299.tiff"),width = 4,height = 2,device = "tiff")

##### combine invasion of A549 and H1299 into one ----
H1299_invasion_summary %>%
  rbind(A549_invasion_summary) %>%
  ggplot(aes(x=group, y=invasion, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=invasion-sd, ymax=invasion+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ x, strip.position = "left",ncol=1) +
  geom_text(aes(y = invasion+sd+5, label=p_labe), size = 5, angle =90) +
  # ylim(0,1.5)+
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  scale_y_continuous(position = "right") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank() ,
    axis.ticks.y = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.x = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text.y = element_text(size=12, angle = 90),
    # legend.key.height = unit(0.15,"inches"),
    # legend.key.width = unit(0.15,"inches"),
    legend.background = element_blank()
  )  +
  ylab(paste("Cells invasion", "(% of Control)")) +
  # labs(title = "Invasion") +
  coord_flip()
ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299-A549.pdf"),width = 4,height = 4,device = "pdf")
ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299-A549.tiff"),width = 4,height = 4,device = "tiff")

ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299-A549-flip.tiff"),width = 4,height = 4,device = "tiff")
##################################################
### CBX2 and EZH2 siRNA cell variability in A549
CBX2_viability_bar_A549 <- readr::read_tsv(file.path(data_path,"2.MTT_A549_value.txt")) %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="viability") 
CBX2_viability_SD_A549 <- readr::read_tsv(file.path(data_path,"2.MTT_A549_STD.txt"))  %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="sd") 
CBX2_viability_P_A549 <- readr::read_tsv(file.path(data_path,"2.MTT_A549_P.txt"))  %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="p") 

CBX2_viability_bar_A549 %>%
  dplyr::inner_join(CBX2_viability_SD_A549,by=c("siRNA","Time (hours)")) %>%
  dplyr::left_join(CBX2_viability_P_A549,by=c("siRNA","Time (hours)")) %>%
  dplyr::mutate(p=as.numeric(p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> CBX2_viability_summary_A549

# ggplot
## bar plot
CBX2_viability_summary_A549<-within(CBX2_viability_summary_A549,
                                    siRNA <- factor(siRNA,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(CBX2_viability_summary_A549,levels(siRNA))
ggplot(CBX2_viability_summary_A549, aes(x=siRNA, y=viability, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~`Time (hours)`,strip.position = "bottom",nrow=1) +
  geom_text(aes(y=viability+sd+0.05,label=p_labe),size=5) +theme_classic() +
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
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot-A549.pdf"),width = 7,height = 3,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot-A549.tiff"),width = 7,height = 3,device = "tiff")

## broken line
CBX2_viability_summary_A549 %>%
  dplyr::mutate(`Time (hours)`=as.character(`Time (hours)`)) %>%
  ggplot(aes(x=`Time (hours)`, y=viability, color=siRNA)) +
  geom_line(aes(group=siRNA)) +
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

ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline-A549.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline-A549.pdf"),width = 4,height = 2,device = "pdf")


##################################################
### CBX2 and EZH2 siRNA cell viability in H1299
CBX2_viability_bar_H1299 <- readr::read_tsv(file.path(data_path,"3.MTT_H1299_value.txt")) %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="viability") 
CBX2_viability_SD_H1299 <- readr::read_tsv(file.path(data_path,"3.MTT_H1299_STD.txt"))  %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="sd") 
CBX2_viability_P_H1299 <- readr::read_tsv(file.path(data_path,"3.MTT_H1299_P.txt"))  %>%
  tidyr::gather(-`Time (hours)`,key = "siRNA",value="p") 

CBX2_viability_bar_H1299 %>%
  dplyr::inner_join(CBX2_viability_SD_H1299,by=c("siRNA","Time (hours)")) %>%
  dplyr::left_join(CBX2_viability_P_H1299,by=c("siRNA","Time (hours)")) %>%
  dplyr::mutate(p=as.numeric(p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> CBX2_viability_summary_H1299

# ggplot
## bar plot
CBX2_viability_summary_H1299<-within(CBX2_viability_summary_H1299,
                                    siRNA <- factor(siRNA,levels = c("Control","siCBX2","siEZH2","siCBX2+EZH2")))
with(CBX2_viability_summary_H1299,levels(siRNA))
ggplot(CBX2_viability_summary_H1299, aes(x=siRNA, y=viability, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~`Time (hours)`,strip.position = "bottom",nrow=1) +
  geom_text(aes(y=viability+sd+0.05,label=p_labe),size=5) +theme_classic() +
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
  xlab("Time (hours)") +
  labs(title = "H1299 cell viability")
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot-H1299.pdf"),width = 7,height = 3,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot-H1299.tiff"),width = 7,height = 3,device = "tiff")

## broken line
CBX2_viability_summary_H1299 %>%
  dplyr::mutate(`Time (hours)`=as.character(`Time (hours)`)) %>%
  ggplot(aes(x=`Time (hours)`, y=viability, color=siRNA)) +
  geom_line(aes(group=siRNA)) +
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

ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline-H1299.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline-H1299.pdf"),width = 4,height = 2,device = "pdf")

##################################################
### CBX2 and EZH2 siRNA cell apoptosis
CBX2_apoptosis_bar.A549 <- readr::read_tsv(file.path(data_path,"5.apoptosis_A549.txt")) %>%
  tidyr::gather(key="group",value="Apoptosis") %>%
  dplyr::mutate(x="A549 cell")
CBX2_apoptosis_bar.H1299 <- readr::read_tsv(file.path(data_path,"5.apoptosis_H1299.txt")) %>%
  tidyr::gather(key="group",value="Apoptosis") %>%
  dplyr::mutate(x="H1299 cell")

data_summary(CBX2_apoptosis_bar.A549,varname = "Apoptosis",groupnames = "group") -> apoptosis_A549_summary
data_summary(CBX2_apoptosis_bar.H1299,varname = "Apoptosis",groupnames = "group") -> apoptosis_H1299_summary

CBX2_apoptosis_bar.A549 %>%
  dplyr::mutate(x=c(rep(c(1,2,3),4))) %>%
  tidyr::spread(key="group",value="Apoptosis") %>%
  as.matrix() -> apoptosis_bar.A549.for_ttest
apoptosis_A549_p <- data.frame(group=apoptosis_A549_summary$group,p=NA)
for(i in 3:5){
  t.test(apoptosis_bar.A549.for_ttest[,2],apoptosis_bar.A549.for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> apoptosis_A549_p[i-1,2]
}

CBX2_apoptosis_bar.H1299 %>%
  dplyr::mutate(x=c(rep(c(1,2,3),4))) %>%
  tidyr::spread(key="group",value="Apoptosis") %>%
  as.matrix() -> apoptosis_bar.H1299.for_ttest
apoptosis_H1299_p <- data.frame(group=apoptosis_H1299_summary$group,p=NA)
for(i in 3:5){
  t.test(apoptosis_bar.H1299.for_ttest[,2],apoptosis_bar.H1299.for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> apoptosis_H1299_p[i-1,2]
}

apoptosis_A549_summary %>%
  dplyr::inner_join(apoptosis_A549_p,by="group") %>%
  dplyr::mutate(x="A549 cell") %>%
  rbind(
    apoptosis_H1299_summary %>%
      dplyr::inner_join(apoptosis_H1299_p,by="group") %>%
      dplyr::mutate(x="H1299 cell")
  ) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) ->apoptosis_combine

# ggplot
## bar plot
apoptosis_combine<-within(apoptosis_combine,
                          group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")))
with(apoptosis_combine,levels(group))

ggplot(apoptosis_combine, aes(x=group, y=Apoptosis, fill = group)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=Apoptosis-sd, ymax=Apoptosis+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~x,strip.position = "left",ncol =1) +
  geom_text(aes(y=Apoptosis+sd+1,label=p_labe),size=5, angle = 90) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text.y = element_text(size = 12,colour = "black", angle = 90)
  ) +
  ylab("Apoptotic cells (%)") +
  # labs(title = "Apoptosis") +
  coord_flip()
ggsave(file.path(result_path,"apoptosis_siCBX-EZH2-H1299-A549.pdf"),width = 6,height = 2,device = "pdf")
ggsave(file.path(result_path,"apoptosis_siCBX-EZH2-H1299-A549.tiff"),width = 6,height = 2,device = "tiff")

ggsave(file.path(result_path,"apoptosis_siCBX-EZH2-H1299-A549-filp.tiff"),width = 4,height = 4,device = "tiff")

##################################################
### siPPAR qPCR
PPAR_bar <- readr::read_tsv(file.path(data_path,"6.PPARg-qpcr.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<0.051 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.051 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(key = "qPCR") 
# ggplot
## bar plot
ggplot(PPAR_bar, aes(x=siRNA, y=mean, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~key, strip.position = "bottom") +
  geom_text(aes(y=mean+SD+0.05,label=p_labe),size = 5) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 10),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size=12),
    # legend.key.height = unit(0.1,"inches"),
    # legend.key.width = unit(0.1,"inches"),
    legend.background = element_blank()
  ) +
  ylab(paste("Relative expression",sep="\n"))
ggsave(file.path(result_path,"PPAR_qPCR.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"PPAR_qPCR.tiff"),width = 3,height = 2,device = "tiff")

