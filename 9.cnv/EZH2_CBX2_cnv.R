# E zhou path 
data_path<- "H:/data/GSCALite/TCGA/cnv"
result_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure"

# HUST path
data_path<- "G:/data/GSCALite/TCGA/cnv"
result_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure"

library(magrittr)
library(ggplot2)


luad_cnv <- readr::read_rds(file.path(data_path,"pancan34_cnv_percent.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()


genelist <- c("EZH2","CBX2")
genelist <- c("CBX6","CBX7")

# filter ------------------------------------------------------------------

luad_cnv %>%
  dplyr::filter(symbol %in% genelist) %>%
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
    legend.position = c(0.8,0.8),
    legend.background = element_blank(),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 20, colour = "black")
  ) ->p;p
ggsave(file.path(result_path,"Figure2","EZH2_CBX2_CNV_percent.pdf"),device = "pdf",width = 5,height = 3)
ggsave(file.path(result_path,"Figure2","EZH2_CBX2_CNV_percent.tiff"),device = "tiff",width = 5,height = 3)

