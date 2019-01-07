library(magrittr)
library(ggplot2)
# data path ---------------------------------------------------------------
# E Zhou -----
data_path <- "H:/data/TCGA/LUAD_methylation"

luad_meth <- readr::read_tsv(file.path(data_path,"LUAD_paired_methylation.txt"))

gene_list=c("CBX2","CBX4","CBX6","CBX7","CBX8")

cbx2 <- readr::read_tsv(file.path("H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/图","CBX2_promoter.txt")) %>%
  dplyr::select(tag)

cbx4 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "CBX4") %>%
  dplyr::arrange(desc(Genomic_Coordinate)) %>% 
  dplyr::filter(Genomic_Coordinate>=Genomic_Coordinate[8]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[3]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
cbx8 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "CBX8") %>%
  dplyr::arrange(desc(Genomic_Coordinate)) %>% 
  dplyr::filter(Genomic_Coordinate>=Genomic_Coordinate[6]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
cbx7 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "CBX7") %>%
  dplyr::arrange(desc(Genomic_Coordinate)) %>% 
  dplyr::filter(Genomic_Coordinate>=Genomic_Coordinate[10]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[2]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
cbx6 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "CBX6") %>%
  dplyr::arrange(desc(Genomic_Coordinate)) %>% 
  dplyr::filter(Genomic_Coordinate>=Genomic_Coordinate[9]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)

rbind(cbx2,cbx4,cbx6,cbx7,cbx8) -> tags

# data filter -------------------------------------------------------------
fn_rank <- function(data){
  data %>%
    dplyr::arrange(Genomic_Coordinate) %>%
    dplyr::mutate(Genomic_Direction=1:nrow(data)) %>%
    dplyr::mutate(Genomic_Direction=as.character(Genomic_Direction))  -> .out
  return(.out)
}
fn_median <- function(data){
  median(data$methy) -> .out
  return(.out)
}

# get promoter tags from mEXPRESS ------
luad_meth %>%
  dplyr::filter(tag %in% tags$tag) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::mutate(group=ifelse(substr(sample,6,6)==1,"N","T")) -> gene_list.methy
luad_meth %>%
  dplyr::filter(tag %in% tags$tag) %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::select(tag,Genomic_Coordinate,Gene_Symbol) %>%
  dplyr::mutate(Genomic_Coordinate=as.numeric(Genomic_Coordinate)) %>%
  tidyr::nest(-Gene_Symbol) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=purrr::map(data,fn_rank)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=paste("P",Genomic_Direction,sep="_"))-> gene_list.tag_posi
gene_list.methy %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group,-Gene_Symbol) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::ungroup()-> gene_list.median.methy

# calculation -------------------------------------------------------------

# EZH2 and CBX2 -----------------------------------------------------------
gene_list.methy %>%
  tidyr::drop_na(methy) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  dplyr::group_by(tag) %>%
  dplyr::do(
    broom::tidy(
      wilcox.test(methy ~ group, data = .)
    )
  ) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.01,"**",sig)) %>%
  dplyr::inner_join(gene_list.methy,by="tag") %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::mutate(laby=max(as.numeric(methy))+0.05) %>%
  dplyr::select(tag,Gene_Symbol,sig,Genomic_Direction,laby) %>%
  dplyr::ungroup() %>%
  unique()-> gene_list.ttest

gene_list.methy %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::drop_na(methy) %>%
  ggplot(aes(x=Genomic_Direction,y=methy,color=group)) +
  geom_boxplot() +
  # geom_jitter(position = "jitter",size=0.5) +
  geom_text(data=gene_list.ttest,mapping=aes(x=Genomic_Direction,y=laby,label=sig),color="black") +
  scale_color_manual(values = c("blue", "red")) +
  # scale_x_discrete(limit = EZH2_CBX2.tag_posi$Genomic_Direction[1:27]) +
  facet_grid(~Gene_Symbol,scales = "free",space = "free") +
  ylab(latex2exp::TeX("Methylation ($\\beta$ value)")) +
  xlab("Promoter") +
  theme(
    axis.line = element_line(color = "black"),
    panel.background  = element_rect(fill = "white", color = "black"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    # axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20),
    axis.ticks.x = element_blank(),
    # legend.position = c(0.9,0.8),
    strip.background = element_rect(fill = "white", colour = "black")
  ) -> EZH2_CBX2.box;EZH2_CBX2.box
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX4678/"

ggsave(file.path(out_path_fig,"CBX_promoter_methylation.pdf"),device = "pdf",width = 6,height = 4)
ggsave(file.path(out_path_fig,"CBX_promoter_methylation.tiff"),device = "tiff",width = 6,height = 4)

#################################
# CNV analysis ---
#################################

data_path<- "H:/data/GSCALite/TCGA/cnv"
luad_cnv <- readr::read_rds(file.path(data_path,"pancan34_cnv_percent.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()
luad_cnv %>%
  dplyr::filter(symbol %in% gene_list) %>%
  dplyr::select(symbol,a_homo,d_homo) %>%
  tidyr::gather(-symbol,key="type",value="CNV") %>%
  dplyr::mutate(Percent=CNV*100) %>%
  dplyr::mutate(type=ifelse(type=="a_homo","Amplification","Deletion"))-> plot_ready

plot_ready %>%
  ggplot(aes(x=symbol,y=Percent,fill=type)) +
  geom_col(position = "stack", width = 0.5) +
  guides(fill=guide_legend(title = "")) +
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  scale_y_continuous(limits = c(0,5)) +
  ylab("CNV Percent (%)") +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    legend.position = c(0.8,0.7),
    legend.background = element_blank(),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 20, colour = "black")
  ) ->p;p
ggsave(file.path(out_path_fig,"CBX_CNV_percent.pdf"),device = "pdf",width = 5,height = 3)
ggsave(file.path(out_path_fig,"CBX_CNV_percent.tiff"),device = "tiff",width = 5,height = 3)
