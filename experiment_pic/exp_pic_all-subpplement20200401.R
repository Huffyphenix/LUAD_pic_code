library(ggplot2)
library(magrittr)

# data path ---------------------------------------------------------------
# Home
data_path <- "F:/胡斐斐/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2020-04-01补实验/data_process"
result_path <- "F:/胡斐斐/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2020-04-01补实验/Pic_by_R"

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
### PPARG chip array
data_Chip <- readr::read_tsv(file.path(data_path,"pparg_chip.txt")) 

data_Chip %>%
  dplyr::select(antibody,sort) %>%
  unique() -> data_Chip_anti_sort


data_Chip_summary <- data_summary(data_Chip, varname="per_input", 
                                  groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(data_Chip_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = ""))


#  bar plot
# get p signif
data_Chip_p <- readr::read_tsv(file.path(data_path,"pparg_chip_pvalue.txt"))
data_Chip_summary %>%
  dplyr::left_join(data_Chip_p,by=c("antibody","target")) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) -> data_Chip_summary

# ggplot
library(ggplot2)
ggplot(data_Chip_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  geom_text(aes(group = targets, y = per_input+sd+(per_input/20), label = p_labe),size=5) +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("IgG", "CBX2", "H2AK119ub","EZH2","H3K27me3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.5,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #分面字体
  ) +
  ylab("% of input") +
  ylim(0,0.35)

ggsave(file.path(result_path,"PPARG_CHIP.pdf"),width = 2,height = 3,device = "pdf")
ggsave(file.path(result_path,"PPARG_CHIP.tiff"),width = 2,height = 3,device = "tiff")

##################################################
### tranwell_supplement
tranwell_supplement <- readr::read_tsv(file.path(data_path,"tranwell_supplement.txt")) %>%
  tidyr::gather(-class,-X1,key="group",value="value") %>%
  tidyr::spread(key="X1",value="value") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe))

# ggplot
tranwell_supplement %>%
  ggplot(aes(x=group, y=`Percent area wound healed`, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=`Percent area wound healed`-sd, ymax=`Percent area wound healed`+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ class,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  geom_text(aes(group = class, y = `Percent area wound healed`+sd+0.001, label = p_labe),size=5) +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#87CEFA", "#458B74"),
                    labels = c("Control", "siCBX2", "siCBX2+siEZH2","siEZH2")) +
  scale_x_discrete(limits = c("Control", "siCBX2","siEZH2", "siCBX2+siEZH2")) +
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

ggsave(file.path(result_path,"transwell_supplement_scratch_test.pdf"),width = 10,height = 3,device = "pdf")
ggsave(file.path(result_path,"transwell_supplement_scratch_test.tiff"),width = 10,height = 3,device = "tiff")
