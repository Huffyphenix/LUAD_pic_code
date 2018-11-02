library(ggplot2)

# data path ---------------------------------------------------------------
# HUST 
data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/原图/data_process"
result_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/实验图片/原图/Pic_by_R"

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
### Chip for CBX2 and EZH2 target RBMS3 and CLDN11
data_Chip <- readr::read_tsv(file.path(data_path,"1.EZH2-CBX2-Chip-CLDN11-RBMS3.txt")) 
data_Chip %>%
  dplyr::select(antibody,sort) %>%
  unique() -> data_Chip_anti_sort


data_Chip_summary <- data_summary(data_Chip, varname="per_input", 
                    groupnames=c("antibody", "targets")) %>%
  dplyr::mutate(targets = as.factor(targets))%>%
  dplyr::inner_join(data_Chip_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = ""))
  
#  bar plot
# get p signif
data_Chip_p <- readr::read_tsv(file.path(data_path,"1.EZH2-CBX2-Chip-CLDN11-RBMS3-pvalue.txt"))

data_Chip_summary %>%
  dplyr::left_join(data_Chip_p,by=c("antibody","targets")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) -> data_Chip_summary

# ggplot
ggplot(data_Chip_summary, aes(x=antibody, y=per_input, fill=antibody)) +
geom_bar(stat="identity", color="black",
         position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ targets,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  geom_text(aes(group = targets, y = per_input+sd+0.01, label = p_labe),size=5) +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("IgG", "CBX2", "H2AK119ub","EZH2","H3K27me3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.4,0.8),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) +
  ylab("%input")

ggsave(file.path(result_path,"CLDN11_RBMS3_CHIP.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"CLDN11_RBMS3_CHIP.tiff"),width = 4,height = 3,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA EDU
EDU <- readr::read_tsv(file.path(data_path,"2.CBX2_EZH2_EDU.txt"),col_names = F) %>%
  tidyr::gather(-X1,key = "group",value="Relative_mRNA_level") %>%
  dplyr::select(-group) %>%
  dplyr::rename("group" = "X1")

data_summary(EDU,varname = "Relative_mRNA_level",groupnames = "group") -> EDU_summary

#  bar plot
# get p signif
EDU %>%
  dplyr::mutate(x=c(rep(1,5),rep(2,5),rep(3,5))) %>%
  tidyr::spread(key="group",value="Relative_mRNA_level") %>%
  as.matrix() -> EDU_for_ttest
edu_p <- data.frame(group=EDU_summary$group,p=NA)
for(i in 3:6){
  t.test(EDU_for_ttest[,2],EDU_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> edu_p[i-1,2]
}

EDU_summary %>%
  dplyr::inner_join(edu_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) -> EDU_summary

# ggplot
ggplot(EDU_summary, aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y = Relative_mRNA_level+sd+0.5, label=p_labe)) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12)
  ) +
  ylab("EdU positive cells(%)")

ggsave(file.path(result_path,"EDU_siCBX-EZH2.pdf"),width = 4,height = 2,device = "pdf")
ggsave(file.path(result_path,"EDU_siCBX-EZH2.tiff"),width = 4,height = 2,device = "tiff")


##################################################
### CBX2 and EZH2 siRNA invasion
CBX2_invasion <- readr::read_tsv(file.path(data_path,"6. CBX2-EZH2-invasion.txt")) %>%
  tidyr::gather(key = "group",value="invasion")

data_summary(CBX2_invasion,varname = "invasion",groupnames = "group") -> CBX2_invasion_summary

#  bar plot
# get p signif
readr::read_tsv(file.path(data_path,"6. CBX2-EZH2-invasion.txt")) %>% as.matrix() -> CBX2_invasion_for_ttest
CBX2_invasion_p <- data.frame(group=CBX2_invasion_summary$group,p=NA)
for(i in 2:5){
  t.test(CBX2_invasion_for_ttest[,1],CBX2_invasion_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> CBX2_invasion_p[i,2]
}

CBX2_invasion_summary %>%
  dplyr::inner_join(CBX2_invasion_p,by="group") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) -> CBX2_invasion_summary

# ggplot
ggplot(CBX2_invasion_summary, aes(x=group, y=invasion, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin=invasion-sd, ymax=invasion+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y = invasion+sd+0.01, label=p_labe)) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank() ,
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.background = element_blank(),
    # legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 10),
    legend.key.width=unit(0.15,"inches"),  # legend size
    legend.key.height=unit(0.15,"inches")
  ) +
  ylab("Flod of invasion")

ggsave(file.path(result_path,"invasion_siCBX-EZH2.pdf"),width = 3.5,height = 2,device = "pdf")
ggsave(file.path(result_path,"invasion_siCBX-EZH2.tiff"),width = 3.5,height = 2,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA Ecell variability
CBX2_viability_bar <- readr::read_tsv(file.path(data_path,"3. CBX2_EZH2_viability.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="viability") 
CBX2_viability_SD <- readr::read_tsv(file.path(data_path,"3. CBX2_EZH2_viability_SD.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="sd") 
CBX2_viability_P <- readr::read_tsv(file.path(data_path,"3. CBX2_EZH2_viability_p.txt")) %>%
  tidyr::gather(-siRNA,key = "time",value="p")

CBX2_viability_bar %>%
  dplyr::inner_join(CBX2_viability_SD,by=c("siRNA","time")) %>%
  dplyr::left_join(CBX2_viability_P,by=c("siRNA","time")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe)) -> CBX2_viability_summary

# ggplot
## bar plot
ggplot(CBX2_viability_summary, aes(x=time, y=viability, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(label=p_labe)) +theme_classic() +
  scale_fill_manual(values=c("black", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank()
  ) +
  ylab("OD 450 nm")
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot-2.pdf"),width = 7,height = 3,device = "pdf")

## broken line
ggplot(CBX2_viability_summary, aes(x=time, y=viability, color=siRNA)) +
  geom_line(aes(group=siRNA)) +
  geom_point() +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.1) +
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
    axis.ticks = element_blank()
  )+
  ylab("OD 450 nm")

ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline.pdf"),width = 3,height = 2,device = "pdf")

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
  geom_bar(stat="identity", color="black", width = 0.5,
           position=position_dodge())  +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y=value+sd+0.5,label=p_labe)) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4"),
                    labels = c("Control", "CBX2 sgRNA")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = c(0.7,0.6),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  ylab("Relative photon flux (X104)")
ggsave(file.path(result_path,"CBX-matas.pdf"),width = 3,height = 3,device = "pdf")
ggsave(file.path(result_path,"CBX-matas.tiff"),width = 3,height = 2.8,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA Ecell apoptosis
CBX2_apoptosis_bar <- readr::read_tsv(file.path(data_path,"4. CBX2_EZH2_apoptosis.txt")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05 ,"ns",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(is.na(p),"",p_labe))

# ggplot
## bar plot
ggplot(CBX2_apoptosis_bar, aes(x=siRNA, y=Apoptosis, fill = siRNA)) + 
  geom_bar(stat="identity", color="black", width = 0.5,
           position=position_dodge())  +
  geom_errorbar(aes(ymin=Apoptosis-sd, ymax=Apoptosis+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y=Apoptosis+sd+0.5,label=p_labe)) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "right",
    legend.title = element_blank()
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
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
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
  geom_text(aes(group = targets, y = mean+sd+0.1, label = p_labe),size=5) +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#87CEFA", "#8B6914"),
                    label = c("Control","siCBX2","siEZH2","siCBX2+siEZH2")
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.7,0.85),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) +
  ylab("Relative mRNA level")

ggsave(file.path(result_path,"siRNA_CBX2_TSG.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"siRNA_CBX2_TSG.tiff"),width = 4,height = 3,device = "tiff")

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
  geom_text(aes(y=Relative_mRNA_level+sd+0.02,label=p_labe)) +theme_classic() +
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
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab("Relative mRNA level")
ggsave(file.path(result_path,"mirna_mimic_barplot.pdf"),width = 4,height = 3,device = "pdf")
ggsave(file.path(result_path,"mirna_mimic_barplot.tiff"),width = 4,height = 3,device = "tiff")

