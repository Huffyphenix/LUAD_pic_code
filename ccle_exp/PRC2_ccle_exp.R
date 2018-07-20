library(magrittr)

# data path ----
ccle_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/ccle/"


# load data ---------------------------------------------------------------

lung_nsc_RPKM <- readr::read_tsv(file.path(ccle_path,"lung_nsc_RPKM"))

genelist <- c("JARID2","AEBP2","EED","SET","SUZ12","RBBP4","RBBP7","EZH1","EZH2") #PRC2 complex



# data fiter --------------------------------------------------------------

lung_nsc_RPKM %>%
  dplyr::filter(Description %in% genelist) -> PRC2_lung_cell_exp


# draw pic ----------------------------------------------------------------
library(ggplot2)

b <- runif(nrow(PRC2_lung_cell_exp), -0.2, 0.2)

PRC2_lung_cell_exp %>%
  ggplot(aes(x=Description,y=exp)) +
  geom_boxplot(fill=NA) +
  geom_point(aes(color=cell_line))
