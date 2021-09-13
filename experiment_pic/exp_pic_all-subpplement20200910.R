.libPaths("J:/library")
library(ggplot2)
library(magrittr)

# data path ---------------------------------------------------------------
data_path <- "F:/胡斐斐/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2020-0804-09补实验/data_process"
result_path <- "F:/胡斐斐/我的坚果云/ENCODE-TCGA-LUAD/实验图片/2020-0804-09补实验/Pic_by_R"

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable to be summariezed
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
  data_sum <- rename(data_sum, c("mean"=varname))
  return(data_sum)
}

##################################################### siCBX2 effect -----
siCBX2_effect <- readr::read_tsv(file.path(data_path,"1.siCBX2-effect.txt")) 
siCBX2_effect %>%
  dplyr::select(antibody,sort) %>%
  unique() -> siCBX2_effect_anti_sort

siCBX2_effect_summary <- data_summary(siCBX2_effect, varname=c("per_input"), 
                                  groupnames=c("antibody")) %>%
  dplyr::inner_join(siCBX2_effect_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::rename("per_input"="varname")

siCBX2_effect %>%
  dplyr::mutate(x=c(1,2,3,1,2,3)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  as.matrix() -> siCBX2_effect_for_ttest
siCBX2_effect_p <- data.frame(antibody=siCBX2_effect_summary$antibody,p=NA)
for(i in 3:ncol(siCBX2_effect_for_ttest)){
  t.test(siCBX2_effect_for_ttest[,2],siCBX2_effect_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> siCBX2_effect_p[i-1,2]
}

siCBX2_effect_summary %>%
  dplyr::inner_join(siCBX2_effect_p,by="antibody") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="A549 cell") %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]}))-> siCBX2_effect_summary

siCBX2_effect_summary <-within(siCBX2_effect_summary,group <- factor(group,levels = c("Control","siCBX2")))
with(siCBX2_effect_summary,levels(group))
ggplot(siCBX2_effect_summary, aes(x=group, y=per_input)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y = per_input+sd+0.1, label=p_labe),size=5) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab("Relative mRNA level")

ggsave(file.path(result_path,"siCBX-effects.pdf"),width = 2,height = 3,device = "pdf")
ggsave(file.path(result_path,"siCBX-effects.tiff"),width = 2,height = 3,device = "tiff")

##################################################### siEZH2 effect ----
siEZH2_effect <- readr::read_tsv(file.path(data_path,"2.siEZH2-effect.txt")) 
siEZH2_effect %>%
  dplyr::select(antibody,sort) %>%
  unique() -> siEZH2_effect_anti_sort

siEZH2_effect_summary <- data_summary(siEZH2_effect, varname=c("per_input"), 
                                      groupnames=c("antibody")) %>%
  dplyr::inner_join(siEZH2_effect_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::rename("per_input"="varname")

siEZH2_effect %>%
  dplyr::mutate(x=c(1,2,3,1,2,3)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  as.matrix() -> siEZH2_effect_for_ttest
siEZH2_effect_p <- data.frame(antibody=siEZH2_effect_summary$antibody,p=NA)
for(i in 3:ncol(siEZH2_effect_for_ttest)){
  t.test(siEZH2_effect_for_ttest[,2],siEZH2_effect_for_ttest[,i]) %>% 
    broom::tidy() %>% .[1,5] -> siEZH2_effect_p[i-1,2]
}

siEZH2_effect_summary %>%
  dplyr::inner_join(siEZH2_effect_p,by="antibody") %>%
  dplyr::mutate(p_labe=ifelse(p<=0.001, "***",p)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.01 & p>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p<=0.05 & p>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p>0.05,"ns",p_labe)) %>%
  dplyr::mutate(x="A549 cell") %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]}))-> siEZH2_effect_summary

siEZH2_effect_summary <-within(siEZH2_effect_summary,group <- factor(group,levels = c("Control","shEZH2")))
with(siEZH2_effect_summary,levels(group))
ggplot(siEZH2_effect_summary, aes(x=group, y=per_input)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(y = per_input+sd+0.1, label=p_labe),size=5) +
  # facet_wrap( ~ targets) +
  # ggpubr::stat_compare_means(ref.group = "Control",method = "t.test",label.y = c(50),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD")) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank() ,
    axis.text = element_text(color = "black",size = 12),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15)
  ) +
  ylab("Relative mRNA level")

ggsave(file.path(result_path,"siEZH2-effects.pdf"),width = 2,height = 3,device = "pdf")
ggsave(file.path(result_path,"siEZH2-effects.tiff"),width = 2,height = 3,device = "tiff")

##################################################### 3.shEZH2-H3K27me3 -----
shEZH2_H3K27me3 <- readr::read_tsv(file.path(data_path,"3.shEZH2-H3K27me3.txt"))

shEZH2_H3K27me3 %>%
  dplyr::select(antibody,sort) %>%
  unique() -> shEZH2_H3K27me3_anti_sort


shEZH2_H3K27me3_summary <- data_summary(shEZH2_H3K27me3, varname="per_input", 
                                  groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(shEZH2_H3K27me3_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]})) 
  
# p value by t test
shEZH2_H3K27me3 %>%
  dplyr::mutate(x=rep(c(1,2,3,4),8)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  tidyr::nest(-target) %>%
  dplyr::mutate(p = purrr::map(data,.f=function(.x){
    broom::tidy(t.test(.x$`Control H3K27me3`,.x$`Control IgG`)) %>% 
      dplyr::mutate(x="Control H3K27me3",y="Control IgG") -> tmp1
    broom::tidy(t.test(.x$`Control H3K27me3`,.x$`shEZH2 H3K27me3`)) %>%
      dplyr::mutate(x="Control H3K27me3",y="shEZH2 H3K27me3") -> tmp2
    broom::tidy(t.test(.x$`shEZH2 H3K27me3`,.x$`shEZH2 IgG`)) %>%
              dplyr::mutate(x="shEZH2 H3K27me3",y="shEZH2 IgG") -> tmp3
    rbind(tmp1,tmp2,tmp3)
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.001, "***",p.value)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.01 & p.value>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.05 & p.value>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value>0.05,"ns",p_labe)) ->shEZH2_H3K27me3_p
shEZH2_H3K27me3_p %>%
  readr::write_tsv(file.path(result_path,"shEZH2_H3K27me3_p.tsv"))

# ggplot
ggplot(shEZH2_H3K27me3_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control IgG", "shEZH2 IgG", "Control H3K27me3","shEZH2 H3K27me3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # legend.key.height = unit(0.15,"inches"),
    # legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #鍒嗛潰瀛椾綋
  ) +
  ylab("% of input") +
  ylim(0,0.35)

ggsave(file.path(result_path,"shEZH2_H3K27me3.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"shEZH2_H3K27me3.tiff"),width = 6,height = 3,device = "tiff")

##################################################### 4.shEZH2-H2AK119ub -----
shEZH2_H2AK119ub <- readr::read_tsv(file.path(data_path,"4.shEZH2-H2AK119ub.txt")) 

shEZH2_H2AK119ub %>%
  dplyr::select(antibody,sort) %>%
  unique() -> shEZH2_H2AK119ub_anti_sort


shEZH2_H2AK119ub_summary <- data_summary(shEZH2_H2AK119ub, varname="per_input", 
                                        groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(shEZH2_H2AK119ub_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]})) 

# p value by t test
shEZH2_H2AK119ub %>%
  dplyr::mutate(x=rep(c(1,2,3,4),8)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  tidyr::nest(-target) %>%
  dplyr::mutate(p = purrr::map(data,.f=function(.x){
    broom::tidy(t.test(.x$`Control H2AK119ub`,.x$`Control IgG`)) %>% 
      dplyr::mutate(x="Control H2AK119ub",y="Control IgG") -> tmp1
    broom::tidy(t.test(.x$`Control H2AK119ub`,.x$`shEZH2 H2AK119ub`)) %>%
      dplyr::mutate(x="Control H2AK119ub",y="shEZH2 H2AK119ub") -> tmp2
    broom::tidy(t.test(.x$`shEZH2 H2AK119ub`,.x$`shEZH2 IgG`)) %>%
      dplyr::mutate(x="shEZH2 H2AK119ub",y="shEZH2 IgG") -> tmp3
    rbind(tmp1,tmp2,tmp3)
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.001, "***",p.value)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.01 & p.value>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.05 & p.value>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value>0.05,"ns",p_labe)) ->shEZH2_H2AK119ub_p
shEZH2_H2AK119ub_p %>%
  readr::write_tsv(file.path(result_path,"shEZH2_H2AK119ub_p.tsv"))

# ggplot
p1<- ggplot(shEZH2_H2AK119ub_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control IgG", "shEZH2 IgG", "Control H2AK119ub","shEZH2 H2AK119ub")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #鍒嗛潰瀛椾綋
  ) +
  ylab("% of input") 
#install.packages("gg.gap")
library(gg.gap)
p1
p2 =gg.gap(plot = p1,
           segments = c(0.2, 0.4),
           tick_width = c(0.2,0.5), 
           rel_heights = c(0.1, 0, 0.2),
           ylim = c(0, 2.5)
)
p2
ggsave(file.path(result_path,"shEZH2_H2AK119ub.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"shEZH2_H2AK119ub.tiff"),width = 6,height = 3,device = "tiff")

##################################################### 5.siCBX2-H3K27me3 -----
siCBX2_H3K27me3 <- readr::read_tsv(file.path(data_path,"5.siCBX2-H3K27me3.txt")) 

siCBX2_H3K27me3 %>%
  dplyr::select(antibody,sort) %>%
  unique() -> siCBX2_H3K27me3_anti_sort


siCBX2_H3K27me3_summary <- data_summary(siCBX2_H3K27me3, varname="per_input", 
                                         groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(siCBX2_H3K27me3_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]})) 

# p value by t test
siCBX2_H3K27me3 %>%
  dplyr::mutate(x=rep(c(1,2,3,4),8)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  tidyr::nest(-target) %>%
  dplyr::mutate(p = purrr::map(data,.f=function(.x){
    broom::tidy(t.test(.x$`Control H3K27me3`,.x$`Control IgG`)) %>% 
      dplyr::mutate(x="Control H3K27me3",y="Control IgG") -> tmp1
    broom::tidy(t.test(.x$`Control H3K27me3`,.x$`siCBX2 H3K27me3`)) %>%
      dplyr::mutate(x="Control H3K27me3",y="siCBX2 H3K27me3") -> tmp2
    broom::tidy(t.test(.x$`siCBX2 H3K27me3`,.x$`siCBX2 IgG`)) %>%
      dplyr::mutate(x="siCBX2 H3K27me3",y="siCBX2 IgG") -> tmp3
    rbind(tmp1,tmp2,tmp3)
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.001, "***",p.value)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.01 & p.value>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.05 & p.value>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value>0.05,"ns",p_labe)) ->siCBX2_H3K27me3_p
siCBX2_H3K27me3_p %>%
  readr::write_tsv(file.path(result_path,"siCBX2_H3K27me3_p.tsv"))

# ggplot
p1 <- ggplot(siCBX2_H3K27me3_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control IgG", "siCBX2 IgG", "Control H3K27me3","siCBX2 H3K27me3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #鍒嗛潰瀛椾綋
  ) +
  ylab("% of input") 
p2 =gg.gap(plot = p1,
           segments = c(0.1, 0.15),
           tick_width = c(0.02,0.05), 
           rel_heights = c(0.1, 0, 0.1),
           ylim = c(0, 0.25)
)
p2
ggsave(file.path(result_path,"siCBX2_H3K27me3.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"siCBX2_H3K27me3.tiff"),width = 6,height = 3,device = "tiff")

##################################################### 6.shEZH2-H2AK119ub -----
siCBX2_H2AK119ub <- readr::read_tsv(file.path(data_path,"6.siCBX2-H2AK119ub.txt")) 

siCBX2_H2AK119ub %>%
  dplyr::select(antibody,sort) %>%
  unique() -> siCBX2_H2AK119ub_anti_sort


siCBX2_H2AK119ub_summary <- data_summary(siCBX2_H2AK119ub, varname="per_input", 
                                         groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(siCBX2_H2AK119ub_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]}))

# p value by t test
siCBX2_H2AK119ub %>%
  dplyr::mutate(x=rep(c(1,2,3,4),8)) %>%
  dplyr::select(-sort) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  tidyr::nest(-target) %>%
  dplyr::mutate(p = purrr::map(data,.f=function(.x){
    broom::tidy(t.test(.x$`Control H2AK119ub`,.x$`Control IgG`)) %>% 
      dplyr::mutate(x="Control H2AK119ub",y="Control IgG") -> tmp1
    broom::tidy(t.test(.x$`Control H2AK119ub`,.x$`siCBX2 H2AK119ub`)) %>%
      dplyr::mutate(x="Control H2AK119ub",y="siCBX2 H2AK119ub") -> tmp2
    broom::tidy(t.test(.x$`siCBX2 H2AK119ub`,.x$`siCBX2 IgG`)) %>%
      dplyr::mutate(x="siCBX2 H2AK119ub",y="siCBX2 IgG") -> tmp3
    rbind(tmp1,tmp2,tmp3)
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.001, "***",p.value)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.01 & p.value>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.05 & p.value>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value>0.05,"ns",p_labe)) ->siCBX2_H2AK119ub_p
siCBX2_H2AK119ub_p %>%
  readr::write_tsv(file.path(result_path,"siCBX2_H2AK119ub_p.tsv"))

# ggplot
p1 <- ggplot(siCBX2_H2AK119ub_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control IgG", "siCBX2 IgG", "Control H2AK119ub","siCBX2 H2AK119ub")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #鍒嗛潰瀛椾綋
  ) +
  ylab("% of input") 
p2 =gg.gap(plot = p1,
           segments = c(0.5, 0.6),
           tick_width = c(0.2,0.5), 
           rel_heights = c(0.1, 0, 0.1),
           ylim = c(0, 2.5)
)
p2
ggsave(file.path(result_path,"siCBX2_H2AK119ub.pdf"),width = 6,height = 3,device = "pdf")
ggsave(file.path(result_path,"siCBX2_H2AK119ub.tiff"),width = 6,height = 3,device = "tiff")

##################################################### 7.siCBX2_siEZH2.txt -----
siCBX2_siEZH2 <- readr::read_tsv(file.path(data_path,"7.siCBX2_siEZH2.txt")) %>%
  tidyr::gather(-X1,key="target",value="per_input") %>%
  dplyr::rename("antibody"="X1")

tibble::tibble(antibody=c("Control","siCBX2","siEZH2","siCBX2+siEZH2"),sort=c(1,2,3,4)) -> sort_index
siCBX2_siEZH2 %>%
  dplyr::inner_join(sort_index,by="antibody") %>%
  dplyr::select(antibody,sort) %>%
  unique() -> siCBX2_siEZH2_anti_sort


siCBX2_siEZH2_summary <- data_summary(siCBX2_siEZH2, varname="per_input", 
                                         groupnames=c("antibody", "target")) %>%
  dplyr::mutate(targets = as.factor(target))%>%
  dplyr::inner_join(siCBX2_siEZH2_anti_sort,by="antibody") %>%
  dplyr::arrange(sort) %>%
  dplyr::mutate(antibody = paste(sort,"_",antibody,sep = "")) %>%
  dplyr::mutate(group = purrr::map(antibody,.f=function(.x){strsplit(.x,"_")[[1]][2]})) %>%
  dplyr::rename("per_input"="varname")

# p value by t test
siCBX2_siEZH2 %>%
  dplyr::mutate(x=rep(c(1,2,3),20)) %>%
  tidyr::spread(key="antibody",value="per_input") %>%
  tidyr::nest(-target) %>%
  dplyr::mutate(p = purrr::map(data,.f=function(.x){
    broom::tidy(t.test(.x$`Control`,.x$`siCBX2`)) %>% 
      dplyr::mutate(x="Control",y="siCBX2") -> tmp1
    broom::tidy(t.test(.x$`Control`,.x$`siEZH2`)) %>% 
      dplyr::mutate(x="Control",y="siEZH2") -> tmp2
    broom::tidy(t.test(.x$`Control`,.x$`siCBX2+siEZH2`)) %>% 
      dplyr::mutate(x="Control",y="siCBX2+siEZH2") -> tmp3
    broom::tidy(t.test(.x$`Control`,.x$`siCBX2+siEZH2`)) %>% 
      dplyr::mutate(x="siCBX2",y="siCBX2+siEZH2") -> tmp4
    broom::tidy(t.test(.x$`Control`,.x$`siCBX2+siEZH2`)) %>% 
      dplyr::mutate(x="siEZH2",y="siCBX2+siEZH2") -> tmp5
    rbind(tmp1,tmp2,tmp3,tmp4,tmp5)
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.001, "***",p.value)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.01 & p.value>0.001,"**",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value<=0.05 & p.value>0.01,"*",p_labe)) %>%
  dplyr::mutate(p_labe=ifelse(p.value>0.05,"ns",p_labe)) ->siCBX2_siEZH2_p
siCBX2_siEZH2_p %>%
  readr::write_tsv(file.path(result_path,"siCBX2_siEZH2_p.tsv"))

# ggplot
ggplot(siCBX2_siEZH2_summary, aes(x=antibody, y=per_input, fill=antibody)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=per_input-sd, ymax=per_input+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap( ~ target,strip.position = "bottom") +
  # ggpubr::stat_compare_means(ref.group = "1_IgG",method = "t.test",label.y = c(0.25),label = "p.signif") +
  theme_classic() +
  scale_fill_manual(values=c("#FFFFFF", "#7FFFD4", "#458B74", "#87CEFA", "#4F94CD"),
                    labels = c("Control", "siCBX2", "siEZH2","siCBX2+siEZH2")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, colour = "black"),
    legend.position = c(0.8,0.89),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    strip.background = element_rect(colour = "white"),
    strip.text = element_text(size = 15, color = "black") #鍒嗛潰瀛椾綋
  ) +
  ylab("% of input") 

ggsave(file.path(result_path,"siCBX2_siEZH2.pdf"),width = 12,height = 6,device = "pdf")
ggsave(file.path(result_path,"siCBX2_siEZH2.tiff"),width = 12,height = 6,device = "tiff")
