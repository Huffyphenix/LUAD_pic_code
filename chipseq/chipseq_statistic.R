######################################################
# For chipseq data visulization
######################################################

# data path ---------------------------------------------------------------

chip_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
library(magrittr)
library(ggplot2)
# laod data ---------------------------------------------------------------

CBX2_targets <- read.table(file.path(chip_path,"CBX2_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb"))
EZH2_targets <- read.table(file.path(chip_path,"H3K27me3_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-fc2"))

# statistic of genes feature
others <- c("3prime_overlapping_ncrna","antisense","IG_V_gene","misc_RNA","processed_transcript","sense_intronic","sense_overlapping","TEC","TR_V_gene","")
pseudogene <- c("IG_V_pseudogene","polymorphic_pseudogene","processed_pseudogene","TR_V_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene")
CBX2_targets %>%
  dplyr::select(V1:V7) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  .$V7 %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::rename("Type"=".") %>%
  dplyr::mutate(Content="Content") %>%
  dplyr::mutate(Class=ifelse(length(grep("pseudogene",Type))==1,"pseudogene","Others")) %>%
  dplyr::mutate(Class=ifelse(length(grep("protein_coding",Type))==1,"protein_coding",Class)) %>%
  dplyr::mutate(Class=ifelse(length(grep("lincRNA",Type))==1,"lincRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(length(grep("miRNA",Type))==1,"miRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(length(grep("rRNA",Type))==1,"rRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(length(grep("snoRNA",Type))==1,"snoRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(length(grep("snRNA",Type))==1,"snRNA",Class)) -> CBX2_targets_statistic
EZH2_targets %>%
  dplyr::select(V1:V7) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  .$V7 %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename("Type"=".") %>%
  dplyr::mutate(Type=as.character(Type))-> EZH2_targets_statistic

fn_pie <- function(data,y,fill){
  data %>%
    ggplot(aes_string(x='Content',y=y,fill=fill)) +
    geom_bar(stat = 'identity', position = 'stack') +
    coord_polar(theta = 'y')
}
fn_pie(CBX2_targets_statistic,"Freq","Class")
CBX2_targets_statistic %>%
  ggplot(aes(x='Content',y=Freq,color=Type)) +
  coord_polar(theta = 'y')
