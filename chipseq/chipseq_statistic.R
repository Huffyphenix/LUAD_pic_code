######################################################
# For chipseq data visulization
######################################################
.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
# data path ---------------------------------------------------------------

# chip_path <- "F:/?ҵļ?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
chip_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4"
data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
library(magrittr)
library(ggplot2)
# laod data ---------------------------------------------------------------

CBX2_targets <- read.table(file.path(chip_path,"CBX2_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb"))
EZH2_targets <- read.table(file.path(chip_path,"H3K27me3_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-fc2"))
Animal_TF <- read.table(file.path(chip_path,"AnimalTFDB_tf.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls"))
Animal_TF %>%
  dplyr::filter(! V1 %in% enzyme_lsit$Symbol) -> Animal_TF
TF <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_TF_FC2_cpm_30")) %>%
  dplyr::select(gene_id,log2FC) 
progene <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_ProGene_FC2_cpm_30")) %>%
  dplyr::select(Gene_id,log2FC) %>%
  dplyr::rename("gene_id"="Gene_id")
TF %>%
  rbind(progene) -> all_gene_de_info
# statistic of target genes feature
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
  dplyr::mutate(Class=ifelse(Type %in% pseudogene,"pseudogene","Others")) %>%
  dplyr::mutate(Class=ifelse(Type %in% "protein_coding","protein_coding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "lincRNA","lincRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "miRNA","miRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "rRNA","rRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "snoRNA","snoRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "snRNA","snRNA",Class)) %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> CBX2_targets_statistic
CBX2_targets %>%
  dplyr::filter(V7 == "protein_coding") %>%
  dplyr::mutate(V7 = as.character(V7)) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  dplyr::mutate(Class=ifelse(V5 %in% as.character(Animal_TF$V1),"TF","non-TF")) %>%
  .$Class %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(Class=as.character(.)) %>%
  dplyr::mutate(Content="Content") %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> CBX2_targets_proteincoding_statistic
  

EZH2_targets %>%
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
  dplyr::mutate(Class=ifelse(Type %in% pseudogene,"pseudogene","Others")) %>%
  dplyr::mutate(Class=ifelse(Type %in% "protein_coding","protein_coding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "lincRNA","lincRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "miRNA","miRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "rRNA","rRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "snoRNA","snoRNA",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% "snRNA","snRNA",Class)) %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> EZH2_targets_statistic
EZH2_targets %>%
  dplyr::filter(V7 == "protein_coding") %>%
  dplyr::mutate(V7 = as.character(V7)) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  dplyr::mutate(Class=ifelse(V5 %in% as.character(Animal_TF$V1),"TF","non-TF")) %>%
  .$Class %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(Class=as.character(.)) %>%
  dplyr::mutate(Content="Content") %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> EZH2_targets_proteincoding_statistic

fn_pie <- function(data,y,fill,name,rcolor,n){
  RColorBrewer::brewer.pal(n,rcolor) -> color
  data %>%
    ggplot(aes_string(x='Content',y=y,fill=fill)) +
    geom_bar(stat = 'identity', position = 'stack', width = 0.5) +
    coord_polar(theta = 'y') +
    labs(x='',y='') +
    scale_fill_manual(name=name,values=color) +
    theme(
      panel.background = element_blank(),
      title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}
fn_pie(CBX2_targets_statistic,"Per","Type",name="Class(5313)","Paired",8)
ggsave(file.path(out_path,"CBX2_targets_constitute.pdf"),width = 4,height = 4)
fn_pie(EZH2_targets_statistic,"Per","Type","Class(13576)","Paired",8)
ggsave(file.path(out_path,"EZH2_targets_constitute.pdf"),width = 4,height = 4)
fn_pie(CBX2_targets_proteincoding_statistic,"Per","Type","Class(13576)","Pastel2",2)
ggsave(file.path(out_path,"CBX2_targets_proteincoding_constitute.pdf"),width = 3,height = 3)
fn_pie(EZH2_targets_proteincoding_statistic,"Per","Type","Class(13576)","Pastel2",2)
ggsave(file.path(out_path,"EZH2_targets_proteincoding_constitute.pdf"),width = 3,height = 3)

# overlap of CBX2 and EZH2 targets of protein coding genes
library(VennDiagram)
EZH2_targets %>%
  dplyr::filter(V7 == "protein_coding") %>%
  dplyr::mutate(V7 = as.character(V7)) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  .$V5 -> EZH2_targets_proteincoding
CBX2_targets %>%
  dplyr::filter(V7 == "protein_coding") %>%
  dplyr::mutate(V7 = as.character(V7)) %>%
  dplyr::select(V5,V7) %>%
  unique() %>%
  .$V5 -> CBX2_targets_proteincoding
vennplot <- venn.diagram(list(CBX2=CBX2_targets_proteincoding,EZH2=EZH2_targets_proteincoding),
             # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.svg"),
             filename = NULL,
             # imagetype = "svg", 
             col = RColorBrewer::brewer.pal(9,"RdBu")[c(1,9)],
             cat.col = RColorBrewer::brewer.pal(9,"RdBu")[c(1,9)],
             cat.pos=c(0,-5),
             scaled = TRUE,
             ext.text = TRUE,
             ext.line.lwd = 2,
             ext.dist = -0.15,
             ext.length = 0.9,
             ext.pos = -4,
             inverted = TRUE,
             rotation.degree = 45,
             cex = 2.5,
             cat.cex = 2.5)
grid.draw(vennplot)
ggsave(file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),width = 4,height = 4)
# overlap of CBX2 and EZH2 targets of protein coding genes
intersect(EZH2_targets_proteincoding,CBX2_targets_proteincoding) -> EHZ2_CBX2_common_targets

EHZ2_CBX2_common_targets %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(gene_id=as.character(.))%>%
  dplyr::left_join(all_gene_de_info,by="gene_id") %>%
  dplyr::mutate(Class=ifelse(log2FC<0,"Down","Up")) %>%
  dplyr::mutate(Class=ifelse(is.na(log2FC),"non-DE",Class)) -> EHZ2_CBX2_common_targets.DE_info
EHZ2_CBX2_common_targets.DE_info %>%
  readr::write_tsv(file.path(chip_path,"EHZ2_CBX2_common_targets.DE_info"))

EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::mutate(Sum=n()) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=n()/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",n(),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::mutate(Content="Content") %>%
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> EHZ2_CBX2_common_targets.DE_statistic
fn_pie(EHZ2_CBX2_common_targets.DE_statistic,"Per","Type","Class(13576)","Pastel2",2)
ggsave(file.path(out_path,"EHZ2_CBX2_common_targets.DE_statistic.pdf"),width = 3,height = 3)
