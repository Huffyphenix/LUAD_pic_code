library(ggplot2)
library(magrittr)

# data path ---------------------------------------------------------------
# HUST 
data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/process_data"
result_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2019-7-29补实验/Pic_by_R"

# E Zhou 
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
### CBX2 and EZH2 siRNA EDU in A549 cell
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

EDU_A549_summary %>%
  dplyr::inner_join(edu_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="A549 cell")-> EDU_A549_summary

# ggplot
EDU_A549_summary<-within(EDU_A549_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+EZH2")))
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

### CBX2 and EZH2 siRNA EDU in H1299 cell
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

EDU_H1299_summary %>%
  dplyr::inner_join(edu_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="H1299 cell")-> EDU_H1299_summary

# ggplot
EDU_H1299_summary<-within(EDU_H1299_summary,group <- factor(group,levels = c("Control","siCBX2","siEZH2","siCBX2+EZH2")))
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


##### combine EDU of A549 and H1299 into one
EDU_H1299_summary %>%
  rbind(EDU_A549_summary) %>%
  ggplot(aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ x, strip.position = "bottom") +
  geom_text(aes(y = Relative_mRNA_level+sd+1, label=p_labe), size = 5) +
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
  ylab(paste("EdU positive", "cells(%)",sep = "\n")) +
  labs(title = "Edu")
ggsave(file.path(result_path,"EDU_siCBX-EZH2-H1299-A549.pdf"),width = 6,height = 2,device = "pdf")
ggsave(file.path(result_path,"EDU_siCBX-EZH2-H1299-A549.tiff"),width = 6,height = 2,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA invasion in A549 cell
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

data_summary(A549_invasion,varname = "invasion",groupnames = "group") -> H1299_invasion_summary

#  bar plot
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

##### combine invasion of A549 and H1299 into one
H1299_invasion_summary %>%
  rbind(A549_invasion_summary) %>%
  ggplot(aes(x=group, y=invasion, fill=group)) +
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
  ylab(paste("Cells invasion", "(% of Control)",sep="\n")) +
  labs(title = "Invasion")
ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299-A549.pdf"),width = 6,height = 2,device = "pdf")
ggsave(file.path(result_path,"invasion_siCBX-EZH2-H1299-A549.tiff"),width = 6,height = 2,device = "tiff")

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
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> CBX2_viability_summary_A549

# ggplot
## bar plot
CBX2_viability_summary_A549<-within(CBX2_viability_summary_A549,
                                    siRNA <- factor(siRNA,levels = c("Control","siCBX2","siEZH2","siCBX2+EZH2")))
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
### PPARG siRNA cell variability
PPARG_viability_bar <- readr::read_tsv(file.path(data_path,"siPPARG_MTT.txt")) 
PPARG_viability_bar %>%
  tidyr::gather(-siRNA,key = "time",value="viability") -> PPARG_viability_bar.gather

data_Chip_summary <- data_summary(PPARG_viability_bar.gather, varname="viability", 
                                  groupnames=c("siRNA", "time")) %>%
  dplyr::mutate(time = as.factor(time))


PPARG_p <- data.frame(siRNA=data_Chip_summary$siRNA,time=data_Chip_summary$time,p=NA) %>%
  tidyr::spread(key="time",value="p") 
for(name in colnames(PPARG_viability_bar)[2:4]){
  t.test(as.vector(t(PPARG_viability_bar[1:4,name])),as.vector(t(PPARG_viability_bar[5:8,name]))) %>% 
    broom::tidy() %>% .[1,5] -> PPARG_p[2,name]
  t.test(as.vector(t(PPARG_viability_bar[1:4,name])),as.vector(t(PPARG_viability_bar[9:11,name]))) %>% 
    broom::tidy() %>% .[1,5] -> PPARG_p[3,name]
}
PPARG_p %>%
  tidyr::gather(-siRNA,key="time",value = "p") -> PPARG_p

PPARG_viability_bar.gather %>%
  dplyr::group_by(siRNA,time) %>%
  dplyr::mutate(mean = mean(viability)) %>%
  dplyr::select(siRNA,time,mean) %>%
  unique() %>%
  dplyr::ungroup()  -> PPARG_viability_mean


PPARG_viability_mean %>%
  dplyr::inner_join(PPARG_p,by=c("siRNA","time")) %>%
  dplyr::left_join(data_Chip_summary,by=c("siRNA","time")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> PPARG_viability_summary

# ggplot
## bar plot
library(ggplot2)
ggplot(PPARG_viability_summary, aes(x=siRNA, y=mean, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~time,strip.position = "bottom",nrow=1) +
  geom_text(aes(y=viability+sd+0.01,label=p_labe),size=5) +theme_classic() +
  scale_fill_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend.key.width=unit(0.15,"inches"),  # legend size
    # legend.key.height=unit(0.15,"inches"),
    strip.text = element_text(size = 15),
    strip.background = element_rect(colour = "white")
  ) +
  ylab("OD 490 nm")
ggsave(file.path(result_path,"viability_siPPAR_barplot.pdf"),width = 7,height = 3,device = "pdf")

## broken line
ggplot(PPARG_viability_summary, aes(x=time, y=mean, color=siRNA)) +
  geom_line(aes(group=siRNA)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_color_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                     labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width=unit(0.15,"inches"),  # legend size
    legend.key.height=unit(0.15,"inches"),
    axis.ticks = element_blank()
  )+
  ylab("OD 490 nm")

ggsave(file.path(result_path,"viability_siPPAR_brokenline.pdf"),width = 4,height = 2,device = "pdf")

##################################################
### TF siRNA Ecell variability
TF_viability_bar <- readr::read_tsv(file.path(data_path,"10. TF_viability_data.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="viability") 
TF_viability_SD <- readr::read_tsv(file.path(data_path,"10. TF_viability_sd.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="sd") 
TF_viability_P <- readr::read_tsv(file.path(data_path,"10. TF_viability_p.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="p")

TF_viability_bar %>%
  dplyr::left_join(TF_viability_SD,by=c("siRNA","time")) %>%
  dplyr::left_join(TF_viability_P,by=c("siRNA","time")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> TF_viability_summary

# ggplot
## bar plot
TF_viability_summary%>%
  dplyr::filter(! siRNA %in% c("siE2F3-1","siE2F3-2")) %>%
  ggplot(aes(x=siRNA, y=viability, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~time,strip.position = "bottom",nrow=1) +
  geom_text(aes(y=viability+sd+0.02,label=p_labe),size=5) +theme_classic() +
  scale_fill_manual(values=c("black", c("#858585"),"#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width=unit(0.15,"inches"),  # legend size
    legend.key.height=unit(0.15,"inches"),
    strip.text = element_text(size = 15),
    strip.background = element_rect(colour = "white")
  ) +
  ylab("OD 450 nm") # CCDK8
ggsave(file.path(result_path,"viability_siTF_barplot.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"viability_siTF_barplot.tiff"),width = 6,height = 3,device = "tiff")

## broken line
TF_viability_summary%>%
  dplyr::filter(! siRNA %in% c("siE2F3-1","siE2F3-2")) %>%
  ggplot(aes(x=time, y=viability, color=siRNA)) +
  geom_line(aes(group=siRNA)) +
  geom_point() +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.1) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_color_manual(values=c("black", c("#858585"),"#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    axis.ticks = element_blank()
  )+
  ylab("OD 470 nm")

ggsave(file.path(result_path,"viability_siTF_brokenline.pdf"),width = 4,height = 2,device = "pdf")

##################################################
### mirna variability
mirna_mimics <- readr::read_tsv(file.path(data_path,"14.mirna_mimic.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe))

# ggplot
mirna_mimics %>%
  ggplot(aes(x=mimics, y=mean, fill = mimics)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~group,strip.position = "bottom",nrow=1,scales = "free") +
  geom_text(aes(y=mean+sd+mean/50,label=p_labe),size=4) +theme_classic() +
  scale_fill_manual(values=c("black","#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend.key.width=unit(0.15,"inches"),  # legend size
    # legend.key.height=unit(0.15,"inches"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour = "white")
  ) +
  ylab("Relative expresion")
ggsave(file.path(result_path,"mirna_mimic.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"mirna_mimic.tiff"),width = 6,height = 3,device = "tiff")

##################################################
### mirna variability
mirna_viability_bar <- readr::read_tsv(file.path(data_path,"12.mirna_viability.txt")) %>%
  tidyr::gather(-c("mimics","group"),key = "time",value="value") %>%
  tidyr::spread(key="group",value = "value") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe))


# ggplot
## bar plot
mirna_viability_bar %>%
  ggplot(aes(x=mimics, y=mean, fill = mimics)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~time,strip.position = "bottom",nrow=1) +
  geom_text(aes(y=mean+sd+0.05,label=p_labe),size=4) +theme_classic() +
  scale_fill_manual(values=c("black","#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend.key.width=unit(0.15,"inches"),  # legend size
    # legend.key.height=unit(0.15,"inches"),
    strip.text = element_text(size = 15),
    strip.background = element_rect(colour = "white")
  ) +
  ylab("OD 490 nm")
ggsave(file.path(result_path,"viability_mirna_barplot.pdf"),width = 7,height = 3,device = "pdf")
ggsave(file.path(result_path,"viability_mirna_barplot.tiff"),width = 7,height = 3,device = "tiff")

## broken line
mirna_viability_bar %>%
  ggplot(aes(x=time, y=mean, color=mimics)) +
  geom_line(aes(group=mimics)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_color_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    axis.ticks = element_blank()
  )+
  ylab("OD 470 nm")

ggsave(file.path(result_path,"viability_mirna_brokenline.pdf"),width = 4,height = 2,device = "pdf")

##################################################
### mice size
CBX2_mice_size <- readr::read_tsv(file.path(data_path,"14.mice_size.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(rank = c(rep(1,4),rep(2,4))) %>%
  dplyr::mutate(group = paste(rank,group,sep = ""))

# ggplot
## bar plot
ggplot(CBX2_mice_size, aes(x=group, y=average, fill = group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())  +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y=average+sd+10,label=p_labe), size = 5) +
  theme_classic() +
  facet_wrap(~ time, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values=c("black", "#7FFFD4"),
                    labels = c("Control", "CBX2 KO")) +
  theme(
    # axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.15,"inches"),
    legend.key.height = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 12)
  ) +
  ylab(latex2exp::TeX(glue::glue(paste("Tumor volume","(mm^{3})")))) +
  xlab("Time (days)")
ggsave(file.path(result_path,"mice_size.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"mice_size.tiff"),width = 4,height = 2,device = "tiff")

ggplot(CBX2_mice_size, aes(x=time, y=average, color = group)) + 
  geom_line(aes(group=group)) +
  geom_point() +
  # geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.1) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_color_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                     labels = c("Control", "CBX2 KO", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width=unit(0.15,"inches"),  # legend size
    legend.key.height=unit(0.15,"inches"),
    axis.ticks = element_blank()
  ) +
  ylab(latex2exp::TeX(glue::glue(paste("Tumor volume","(mm^{3})")))) +
  xlab("Time (days)")
ggsave(file.path(result_path,"mice_size_broken_line.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"mice_size_broken_line.tiff"),width = 3,height = 2,device = "tiff")

##################################################
### mice weight
mice_weight <- readr::read_tsv(file.path(data_path,"15.mice_weight.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(rank = 1:2) %>%
  dplyr::mutate(group = paste(rank,group,sep = ""))

# ggplot
## bar plot
ggplot(mice_weight, aes(x=group, y=mean, fill = group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y=mean+sd+10,label=p_labe), size = 5) +
  theme_classic() +
  # facet_wrap(~ group, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values=c("black", "#7FFFD4"),
                    labels = c("Control", "CBX2 sgRNA")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.15,"inches"),
    legend.key.height = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 12)
  ) +
  ylab(paste("Tumor weight","(mg)",sep = "\n")) 
ggsave(file.path(result_path,"mice_weight.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"mice_weight.tiff"),width = 3,height = 2,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA Ecell matas
CBX2_mice_matas <- readr::read_tsv(file.path(data_path,"7. mice_matas.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(rank = 1:2) %>%
  dplyr::mutate(siRNA = paste(rank,siRNA,sep = ""))

# ggplot
## bar plot
ggplot(CBX2_mice_matas, aes(x=siRNA, y=value, fill = siRNA)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())  +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y=value+sd+0.5,label=p_labe), size = 5) +
  theme_classic() +
  scale_fill_manual(values=c("black", "#7FFFD4"),
                    labels = c("Control", "CBX2 sgRNA")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.15,"inches"),
    legend.key.height = unit(0.15,"inches")
  ) +
  ylab(paste("Relative photon", "flux (X104)", sep="\n"))
ggsave(file.path(result_path,"CBX-matas.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"CBX-matas.tiff"),width = 3,height = 2,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA Ecell apoptosis
CBX2_apoptosis_bar <- readr::read_tsv(file.path(data_path,"4. CBX2_EZH2_apoptosis.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(x="Apoptosis")

# ggplot
## bar plot
ggplot(CBX2_apoptosis_bar, aes(x=siRNA, y=Apoptosis, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=Apoptosis-sd, ymax=Apoptosis+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~x,strip.position = "bottom") +
  geom_text(aes(y=Apoptosis+sd+1,label=p_labe),size=5) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 12,colour = "black")
  ) +
  ylab("Apoptotic cells (%)")
ggsave(file.path(result_path,"apoptosis_siCBX-EZH2_barplot.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"apoptosis_siCBX-EZH2_barplot.tiff"),width = 3,height = 2,device = "tiff")

##################################################
### siRNA target PPARG
siRNA_target_PPARG <- readr::read_tsv(file.path(data_path,"8. siCBX_PPARG.txt")) 
siRNA_target_PPARG <- readr::read_tsv(file.path(data_path,"8. siEZH_PPARG.txt")) 

siRNA_target_PPARG %>%
  tidyr::gather(-siRNA,key="targets",value="Relative_mRNA_level") -> siRNA_target_PPARG.gather

siRNA_target_PPARG_summary <- data_summary(siRNA_target_PPARG.gather, varname="Relative_mRNA_level", 
                                           groupnames=c("siRNA", "targets")) 

#  bar plot
# get p signif
siRNA_target_PPARG

PPARG_p <- data.frame(siRNA=siRNA_target_PPARG_summary$siRNA,targets=siRNA_target_PPARG_summary$targets,p=NA) %>%
  tidyr::spread(key="targets",value="p") 
for(name in colnames(siRNA_target_PPARG)[2:5]){
  t.test(as.vector(t(siRNA_target_PPARG[1:3,name])),as.vector(t(siRNA_target_PPARG[4:6,name]))) %>% 
    broom::tidy() %>% .[1,5] -> PPARG_p[2,name]
  t.test(as.vector(t(siRNA_target_PPARG[1:3,name])),as.vector(t(siRNA_target_PPARG[7:9,name]))) %>% 
    broom::tidy() %>% .[1,5] -> PPARG_p[3,name]
}
PPARG_p %>%
  tidyr::gather(-siRNA,key="targets",value="p") %>%
  dplyr::left_join(siRNA_target_PPARG_summary,by=c("siRNA","targets")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) -> siRNA_target_PPARG_summary

# ggplot
ggplot(siRNA_target_PPARG_summary, aes(x=siRNA, y=Relative_mRNA_level, fill=siRNA)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ targets,strip.position = "bottom",nrow = 1) +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  geom_text(aes(group = targets, y = Relative_mRNA_level+sd+0.05, label = p_labe),size=5) +
  scale_fill_manual(#values=c("#FFFFFF", "#7FFFD4", "#458B74") # For CBX2
    values=c("#FFFFFF", "#87CEFA", "#4F94CD") # For EZH2
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.2,0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) +
  ylab("Relative mRNA level")

ggsave(file.path(result_path,"siRNA_CBX2_PPAR.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"siRNA_CBX2_PPAR.tiff"),width = 4,height = 3,device = "tiff")
ggsave(file.path(result_path,"siRNA_EZH2_PPAR.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"siRNA_EZH2_PPAR.tiff"),width = 4,height = 3,device = "tiff")


##################################################
### siRNA target TSG
siRNA_target_TSG <- readr::read_tsv(file.path(data_path,"8. sCBX-EZH-TSG.txt")) %>%
  tidyr::gather(-c(targets,group),key="siRNA",value="value") %>%
  tidyr::spread(key="group",value="value") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.053 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.053,"ns",p_labe)) %>%
  dplyr::mutate(sort=rep(c(1,2,4,3),2)) %>%
  dplyr::mutate(siRNA = paste(sort,siRNA,sep = ""))

# ggplot
ggplot(siRNA_target_TSG, aes(x=siRNA, y=mean, fill=siRNA)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ targets,strip.position = "bottom",nrow = 1) +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  geom_text(aes(group = targets, y = mean+sd+0.2, label = p_labe),size=5) +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#87CEFA", "#4F94CD"),
                    label = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")
  ) +
  ylim(0,7)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.4,0.9),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) +
  ylab("Relative mRNA level")

ggsave(file.path(result_path,"siRNA_CBX2_TSG.pdf"),width = 2.5,height = 3,device = "pdf")
ggsave(file.path(result_path,"siRNA_CBX2_TSG.tiff"),width = 2.5,height = 3,device = "tiff")

##################################################
### CBX2 and EZH2 miRNA mimics
miRNA_mimics_bar <- readr::read_tsv(file.path(data_path,"5. miRNA mimics.txt")) %>%
  tidyr::gather(-c(siRNA,group),key="key",value="value") %>%
  tidyr::spread(key="group",value="value") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<0.051 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.051 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe))

# ggplot
## bar plot
ggplot(miRNA_mimics_bar, aes(x=key, y=Relative_mRNA_level, fill = key)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~siRNA, strip.position = "bottom") +
  geom_text(aes(y=Relative_mRNA_level+sd+0.02,label=p_labe),size = 5) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.85),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab("Relative mRNA level")
ggsave(file.path(result_path,"mirna_mimic_barplot.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"mirna_mimic_barplot.tiff"),width = 4,height = 3,device = "tiff")
##################################################
### siPPAR transwell
PPAR_bar <- readr::read_tsv(file.path(data_path,"siPPAR-transwell.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<0.051 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.051 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(key = "Invasion") 
# ggplot
## bar plot
ggplot(PPAR_bar, aes(x=siRNA, y=mean, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~key, strip.position = "bottom") +
  geom_text(aes(y=mean+SD+2,label=p_labe),size = 5) +theme_classic() +
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
  ylab(paste("Cells invasion", "(% of Control)",sep="\n"))
ggsave(file.path(result_path,"PPAR_transwell.pdf"),width = 3,height = 2,device = "pdf")
ggsave(file.path(result_path,"PPAR_transwell.tiff"),width = 3,height = 2,device = "tiff")
##################################################
### siTF transwell
TF_EZH2_bar <- readr::read_tsv(file.path(data_path,"9. siTF-transwell.txt")) %>%
  tidyr::gather(-c(group),key="siRNA",value="value") %>%
  tidyr::spread(key="group",value="value") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<0.051 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.051 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) %>%
  dplyr::mutate(key = "Invasion") %>%
  dplyr::mutate(mean = mean * 100) %>%
  dplyr::mutate(sd = sd * 100)

# ggplot
## bar plot
ggplot(TF_EZH2_bar, aes(x=siRNA, y=mean, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~key, strip.position = "bottom") +
  geom_text(aes(y=mean+sd+2,label=p_labe),size = 5) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = c(0.5,0.9),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 8),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size=12),
    legend.key.height = unit(0.1,"inches"),
    legend.key.width = unit(0.1,"inches"),
    legend.background = element_blank()
  ) +
  ylab(paste("Cells invasion", "(% of Control)",sep="\n"))
ggsave(file.path(result_path,"TF_transwell.pdf"),width = 2.5,height = 2,device = "pdf")
ggsave(file.path(result_path,"TF_transwell.tiff"),width = 2.5,height = 2,device = "tiff")


##################################################
### siTF CBX2 and EZH2 
siTF_targets_exp <- readr::read_tsv(file.path(data_path,"siTF-targets.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<0.051 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.051 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) 

# ggplot
## bar plot
ggplot(siTF_targets_exp, aes(x=siRNA, y=mean, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~targets, strip.position = "bottom") +
  geom_text(aes(y=mean+sd+0.1,label=p_labe),size = 5) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD","#C7C7C7", "#474747")) +
  ylim(0,2.8) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.8),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab("Relative mRNA level")
ggsave(file.path(result_path,"siTF_targets.pdf"),width = 5,height = 3,device = "pdf")
ggsave(file.path(result_path,"siTF_targets.tiff"),width = 6,height = 3,device = "tiff")
