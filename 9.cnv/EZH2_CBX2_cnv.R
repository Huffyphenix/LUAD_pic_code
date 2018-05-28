data_path<- "H:/data/GSCALite/TCGA/cnv"
result_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure"

library(magrittr)
library(ggplot2)


luad_cnv <- readr::read_rds(file.path(data_path,"pancan34_cnv_percent.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()


genelist <- c("EZH2","CBX2")


# filter ------------------------------------------------------------------

luad_cnv %>%
  dplyr::filter(symbol %in% genelist) %>%
  dplyr::select(symbol,a_homo,d_homo) %>%
  tidyr::gather(-symbol,key="type",value="CNV") %>%
  dplyr::mutate(Percent=CNV*100) %>%
  dplyr::mutate(type=ifelse(type=="a_homo","Amplification","Deletion"))-> plot_ready

plot_ready %>%
  ggplot(aes(x=symbol,y=Percent,fill=type)) +
  geom_col(position = "stack") +
  guides(fill=guide_legend(title = "")) +
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  ylab("Percent (%)") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    legend.position = c(0.8,0.9),
    legend.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size = 12)
  ) ->p;p
ggsave(filename = "EZH2_CBX2_CNV_percent.pdf",path = result_path,device = "pdf",width = 4,height = 4)
ggsave(filename = "EZH2_CBX2_CNV_percent.tiff",path = result_path,device = "tiff",width = 4,height = 4)

