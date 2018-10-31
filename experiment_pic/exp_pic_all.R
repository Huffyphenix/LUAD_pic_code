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
comp_list <- list(c("1_IgG","2_CBX2"),c("1_IgG","3_H2AK119ub"),c("1_IgG","4_EZH2"),c("1_IgG","5_H3K27me3"))
comp_list <- list(c("IgG","CBX2"),c("IgG","H2AK119ub"),c("IgG","EZH2"),c("IgG","H3K27me3"))
data_Chip %>%
  ggpubr::ggboxplot(x = "antibody", y = "per_input",
                    fill  = "antibody", palette = "antibody") +
  facet_wrap( ~ targets) +
  # coord_flip() +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "t.test",label.y = c(0.25,0.27,0.29,0.31),label = "p.signif") +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("IgG", "CBX2", "H2AK119ub","EZH2","H3K27me3")) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("%input")

ggsave()

# ggplot
ggplot(data_Chip_summary, aes(x=targets, y=per_input, fill=antibody)) +
geom_bar(stat="identity", color="black",
         position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("IgG", "CBX2", "H2AK119ub","EZH2","H3K27me3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "top",
    legend.direction  = c("horizontal"),
    legend.title = element_blank()
  ) +
  ylab("%input")

ggsave(file.path(result_path,"CLDN11_RBMS3_CHIP.pdf"),width = 5,height = 3,device = "pdf")
ggsave(file.path(result_path,"CLDN11_RBMS3_CHIP.tiff"),width = 5,height = 3,device = "tiff")

##################################################
### CBX2 and EZH2 siRNA EDU
EDU <- readr::read_tsv(file.path(data_path,"2.CBX2_EZH2_EDU.txt"),col_names = F) %>%
  tidyr::gather(-X1,key = "group",value="Relative_mRNA_level") %>%
  dplyr::select(-group) %>%
  dplyr::rename("group" = "X1")

data_summary(EDU,varname = "Relative_mRNA_level",groupnames = "group") -> EDU_summary

#  bar plot
# get p signif
comp_list <- list(c("Control","siEZH2-1"),c("Control","siEZH2-2"),c("Control","siCBX2-1"),c("Control","siCBX2-2"))
EDU %>%
  ggpubr::ggboxplot(x = "group", y = "Relative_mRNA_level",
                    fill  = "group", palette = "group") +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "t.test",label.y = c(50,51,52,53,54),label = "p.signif") +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("IgG", "CBX2", "H2AK119ub","EZH2","H3K27me3")) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("%input")


# ggplot
ggplot(EDU_summary, aes(x=group, y=Relative_mRNA_level, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin=Relative_mRNA_level-sd, ymax=Relative_mRNA_level+sd), width=.2,
                position=position_dodge(.9)) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "top",
    legend.direction  = c("horizontal"),
    legend.title = element_blank()
  ) +
  ylab("%input")

ggsave(file.path(result_path,"EDU_siCBX-EZH2.pdf"),width = 5,height = 3,device = "pdf")


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
  dplyr::mutate(p=ifelse(p>0.05,"ns","*")) %>%
  dplyr::mutate(p=ifelse(p<=0.01,"**",p)) %>%
  dplyr::mutate(p=ifelse(p<=0.001,"***",p)) %>%
  dplyr::mutate(p=ifelse(p<=0.0001,"****",p)) %>%
  dplyr::mutate(p=ifelse(is.na(p),"",p)) -> CBX2_viability_summary

# ggplot
## bar plot
ggplot(CBX2_viability_summary, aes(x=time, y=viability, fill = siRNA)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())  +
  geom_errorbar(aes(ymin=viability-sd, ymax=viability+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(label=p)) +theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2-1", "siCBX2-2","siEZH2-1","siEZH2-2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = "top",
    legend.direction  = c("horizontal"),
    legend.title = element_blank()
  ) +
  ylab("OD 450 nm")
ggsave(file.path(result_path,"viability_siCBX-EZH2_barplot.pdf"),width = 5,height = 3,device = "pdf")

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
    legend.position = "none",
    legend.title = element_blank(),
    axis.ticks = element_blank()
  )+
  ylab("OD 450 nm")

ggsave(file.path(result_path,"viability_siCBX-EZH2_brokenline.pdf"),width = 2,height = 1,device = "pdf")
