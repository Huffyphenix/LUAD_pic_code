############## overlap between each DE gene set (CBX2, EZH2 and both, high vs. low)
library(VennDiagram)
# path --------------------------------------------------------------------
# HUST
basic_path <- file.path("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD")
# exp_path <- "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
res_path <- file.path(basic_path,"Figure/DE_between_high_low_CBX2_EZH2")


# load data ---------------------------------------------------------------

CBX2_DE <- readr::read_tsv(file.path(res_path,"CBX2","DE_gene_res.tsv")) %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::filter(abs(log2FC_HvsL)>=1)

EZH2_DE <- readr::read_tsv(file.path(res_path,"EZH2","DE_gene_res.tsv")) %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::filter(abs(log2FC_HvsL)>=1)

CBX2_EZH2_DE <- readr::read_tsv(file.path(res_path,"CBX2_and_EZH2","DE_gene_res.tsv")) %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::filter(abs(log2FC_HvsL)>=1)


# overlap -----------------------------------------------------------------
# upregulate genes 
CBX2_DE %>%
  dplyr::filter(log2FC_HvsL>0) %>%
  .$symbol -> CBX2_DE.up
EZH2_DE %>%
  dplyr::filter(log2FC_HvsL>0) %>%
  .$symbol -> EZH2_DE.up
CBX2_EZH2_DE %>%
  dplyr::filter(log2FC_HvsL>0) %>%
  .$symbol -> CBX2_EZH2_DE.up
vennplot <- venn.diagram(list(CBX2_up=CBX2_DE.up,
                              EZH2_up=EZH2_DE.up,
                              CBX2_EZH2_up=CBX2_EZH2_DE.up),
                         # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),
                         filename = NULL,
                         # imagetype = "pdf",
                         col = RColorBrewer::brewer.pal(7,"Set1")[c(7,2,1)],
                         cat.col = RColorBrewer::brewer.pal(7,"Set1")[c(7,2,1)],
                         # cat.pos=c(0,-5),
                         scaled = TRUE,
                         ext.text = TRUE,
                         ext.line.lwd = 2,
                         ext.dist = -0.15,
                         ext.length = 0.9,
                         ext.pos = -4,
                         inverted = TRUE,
                         # rotation.degree = 45,
                         cex = 2.5,
                         cat.cex = 2.5)
grid.draw(vennplot)

# downregulate genes
CBX2_DE %>%
  dplyr::filter(log2FC_HvsL<0) %>%
  .$symbol -> CBX2_DE.down
EZH2_DE %>%
  dplyr::filter(log2FC_HvsL<0) %>%
  .$symbol -> EZH2_DE.down
CBX2_EZH2_DE %>%
  dplyr::filter(log2FC_HvsL<0) %>%
  .$symbol -> CBX2_EZH2_DE.down
vennplot <- venn.diagram(list(CBX2_down=CBX2_DE.down,
                              EZH2_down=EZH2_DE.down,
                              CBX2_EZH2_down=CBX2_EZH2_DE.down),
                         # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),
                         filename = NULL,
                         # imagetype = "pdf",
                         col = RColorBrewer::brewer.pal(7,"Set1")[c(7,2,1)],
                         cat.col = RColorBrewer::brewer.pal(7,"Set1")[c(7,2,1)],
                         # cat.pos=c(0,-5),
                         scaled = TRUE,
                         ext.text = TRUE,
                         ext.line.lwd = 2,
                         ext.dist = -0.15,
                         ext.length = 0.9,
                         ext.pos = -4,
                         inverted = TRUE,
                         # rotation.degree = 45,
                         cex = 2.5,
                         cat.cex = 2.5)
grid.draw(vennplot)

