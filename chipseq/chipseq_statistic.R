######################################################
# For chipseq data visulization
######################################################
.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
# data path ---------------------------------------------------------------

# chip_path <- "F:/?ҵļ?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"

chip_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4"
data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
data_path_3 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

# chip_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
# out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4"
# chip_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"
# out_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4"
# data_path<- "H:/data"
# # data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
# data_path_2 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"


library(magrittr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
# laod data ---------------------------------------------------------------

CBX2_targets <- readr::read_tsv(file.path(chip_path,"180426_promoter_add-8-2kb","CBX2_pva_3-summits_OVERLAP_GRCh38-hg38_TSS-8-2kb.fc2"),col_names = F)
EZH2_targets <- readr::read_tsv(file.path(chip_path,"180426_promoter_add-8-2kb","H3K27me3_pva_3-summits_OVERLAP_GRCh38-hg38-TSS-8-2kb"),col_names = F)
# CBX2_targets$X7 %>% unique() %>%
#   bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db) -> CBX2_targets_entrez
# EZH2_targets$X7 %>% unique() %>%
#   bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db) -> EZH2_targets_entrez
# co_occurence_targets_ignore_overlap <- readr::read_tsv(file.path(chip_path,"180426_promoter_add-8-2kb","co_occurence_targets"),col_names = F) %>%
#   .$X1 %>%
#   bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
EZH2_targets_inoverlappeaks <- readr::read_tsv(file.path(chip_path,"180413_peak_overlap_CBX2_H3K27me3","EZH2_summits_overlap_CBX2-OVERLAP-GRCh38-TSS-8-2kb"),col_names = F)
CBX2_targets_inoverlappeaks <- readr::read_tsv(file.path(chip_path,"180413_peak_overlap_CBX2_H3K27me3","CBX2_summits_overlap_H3K27me3-OVERLAP-GRCh38-TSS-8-2kb"),col_names = F)

Animal_TF <-  readr::read_tsv(file.path(data_path,"AnimalTFDB","Homo_sapiens_transcription_factors_gene_list.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls")) %>%
  .$Symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Animal_TF %>%
  dplyr::filter(! Entrez_ID %in% enzyme_lsit$ENTREZID) -> Animal_TF
Animal_TF$Entrez_ID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db) -> Animal_TF_entrez

TF <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_TF_FC2_cpm_30")) %>%
  dplyr::select(gene_id,log2FC) 
progene <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_ProGene_FC2_cpm_30")) %>%
  dplyr::select(Gene_id,log2FC) %>%
  dplyr::rename("gene_id"="Gene_id")
TF %>%
  rbind(progene) -> all_gene_de_info

TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil
all_gene_nofil %>%
  dplyr::filter(log2FC!="NA") %>%
  dplyr::filter(abs(log2FC)>=0.585) %>%
  dplyr::filter(prob>0.90) %>%
  dplyr::filter(case_mean>=30|con_mean>=30) -> all_gene_prob0.99_FC1.5_exp30

# statistic of target genes feature
# biotypes group comes from：ensembl.
# https://asia.ensembl.org/Help/Faq?id=468
# http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
protein_coding <- c("protein_coding","IG_V_gene","TR_V_gene","polymorphic_pseudogene")
Pseudogene <- c("IG_V_pseudogene","TR_V_pseudogene","pseudogene","transcribed_unprocessed_pseudogene","transcribed_processed_pseudogene",
                "unprocessed_pseudogene","processed_pseudogene")
Long_noncoding <- c("3prime_overlapping_ncrna","antisense","lincRNA","processed_transcript","sense_overlapping","sense_intronic")
Short_noncoding <- c("miRNA","misc_RNA","Mt_tRNA","Mt_rRNA","rRNA","snoRNA","snRNA")
# others <- c("3prime_overlapping_ncrna","antisense","IG_V_gene","misc_RNA","processed_transcript","sense_intronic","sense_overlapping","TEC","TR_V_gene","")
# pseudogene <- c("IG_V_pseudogene","polymorphic_pseudogene","processed_pseudogene","TR_V_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene")
CBX2_targets %>%
  # dplyr::select(X1:X7) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X7 %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::rename("Type"=".") %>%
  dplyr::mutate(Content="Content") %>%
  dplyr::mutate(Class=ifelse(Type %in% pseudogene,"pseudogene","Others")) %>%
  dplyr::mutate(Class=ifelse(Type %in% protein_coding,"Protein_coding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% Long_noncoding,"Long_noncoding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% Short_noncoding,"Short_noncoding",Class)) %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> CBX2_targets_statistic
CBX2_targets %>%
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  dplyr::mutate(Class=ifelse(X8 %in% as.character(Animal_TF_entrez$SYMBOL),"TF","non-TF")) %>%
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
  # dplyr::select(V1:V7) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X7 %>%
  table() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::rename("Type"=".") %>%
  dplyr::mutate(Content="Content") %>%
  dplyr::mutate(Class=ifelse(Type %in% pseudogene,"pseudogene","Others")) %>%
  dplyr::mutate(Class=ifelse(Type %in% protein_coding,"Protein_coding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% Long_noncoding,"Long_noncoding",Class)) %>%
  dplyr::mutate(Class=ifelse(Type %in% Short_noncoding,"Short_noncoding",Class)) %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(Per=sum(Freq)/Sum) %>%
  dplyr::mutate(Type=paste(Class,"(",sum(Freq),", ",round(Per*100,1),"%)",sep="")) %>% 
  dplyr::select(Content,Class,Per,Type,Sum) %>%
  unique() -> EZH2_targets_statistic
EZH2_targets %>%
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  dplyr::mutate(Class=ifelse(X8 %in% as.character(Animal_TF_entrez$SYMBOL),"TF","non-TF")) %>%
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
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X8 -> EZH2_targets_proteincoding

EZH2_targets_inoverlappeaks %>%
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X8  -> EZH2_targets_proteincoding_inoverlappeak

CBX2_targets %>%
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X8 -> CBX2_targets_proteincoding
CBX2_targets_inoverlappeaks %>%
  dplyr::filter(X7 %in% protein_coding) %>%
  dplyr::mutate(X7 = as.character(X7)) %>%
  dplyr::select(X7,X8) %>%
  unique() %>%
  .$X8 -> CBX2_targets_proteincoding_inoverlappeak

library(VennDiagram)
vennplot <- venn.diagram(list(CBX2_All=CBX2_targets_proteincoding,
                              EZH2_All=EZH2_targets_proteincoding,
                              CBX2_Co=CBX2_targets_proteincoding_inoverlappeak,
                              EZH2_Co=EZH2_targets_proteincoding_inoverlappeak),
             # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),
             filename = NULL,
             # imagetype = "pdf", 
             col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2,5,1)],
             cat.col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2,5,1)],
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
#ggsave(file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),width = 4,height = 4)

### overlap of CBX2 and EZH2 targets of protein coding genes
intersect(EZH2_targets_proteincoding,CBX2_targets_proteincoding) -> EZH2_CBX2_common_targets
intersect(EZH2_targets_proteincoding_inoverlappeak,CBX2_targets_proteincoding_inoverlappeak) -> EZH2_CBX2_common_peak_common_targets
setdiff(EZH2_targets_proteincoding_inoverlappeak,CBX2_targets_proteincoding_inoverlappeak)


EZH2data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil


# DE info in LUAD of common targets ---------------------------------------


tcga_geneid <- readr::read_tsv("F:/我的坚果云/ENCODE-TCGA-LUAD/TCGA_gene_info/TCGA_all_gene_id.txt")
all_gene_nofil %>%
  dplyr::inner_join(tcga_geneid,by="gene_id") %>%
  dplyr::mutate(entrez_id=as.character(entrez_id)) -> all_gene_nofil.entrez

EZH2_CBX2_common_targets %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db) -> EZH2_CBX2_common_targets.entrez
EZH2_CBX2_common_targets %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::filter(! . %in% EZH2_CBX2_common_targets.entrez$SYMBOL) %>%
  .$. %>%
  bitr(fromType = "ALIAS",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = F) %>%
  dplyr::rename("SYMBOL"="ALIAS") -> EZH2_CBX2_common_targets.nomap.entrez
EZH2_CBX2_common_targets.entrez %>%
  rbind(EZH2_CBX2_common_targets.nomap.entrez) %>%
  dplyr::rename("gene_id"="SYMBOL","entrez_id"="ENTREZID") %>%
  dplyr::mutate(entrez_id=as.character(entrez_id))-> EZH2_CBX2_common_targets.entrez.all
  
EZH2_CBX2_common_targets.entrez.all %>%
  dplyr::as.tbl() %>%
  dplyr::left_join(all_gene_nofil.entrez,by="entrez_id") %>%
  dplyr::mutate(Class=ifelse(log2FC<= (-0.585) & prob>=0.99,"Down","non-DE")) %>%
  dplyr::mutate(Class=ifelse(log2FC>= 0.585 & prob>=0.99,"Up",Class)) %>%
  dplyr::mutate(Class=ifelse(is.na(log2FC),"non-DE",Class)) -> EHZ2_CBX2_common_targets.DE_info
EHZ2_CBX2_common_targets.DE_info %>%
  readr::write_tsv(file.path(chip_path,"common-targets-180426-new","all_EHZ2_CBX2_common_targets.DE_info"))
EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::filter(gene_id.x %in% EZH2_CBX2_common_peak_common_targets) %>%
  readr::write_tsv(file.path(chip_path,"common-targets-180426-new","only_EHZ2_CBX2_common_peaks-targets.DE_info"))
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


# pathway enrichemnt ------------------------------------------------------
library(clusterProfiler)
library(enrichplot)

##  Go enrichment ----
EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::filter(! Class=="non-DE") -> EHZ2_CBX2_common_targets.DE
EHZ2_CBX2_common_targets.DE$entrez_id %>%
  enrichGO(OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE) -> DE_Go_enrichment
heatplot(DE_Go_enrichment, foldChange=EHZ2_CBX2_common_targets.DE$log2FC)

EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::filter(Class=="Down") -> EHZ2_CBX2_common_targets.Down
EHZ2_CBX2_common_targets.Down$entrez_id %>%
  enrichGO(OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE) -> Down_Go_enrichment

EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::filter(Class=="Up") -> EHZ2_CBX2_common_targets.Up
EHZ2_CBX2_common_targets.Down$entrez_id %>%
  enrichGO(OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE) -> Up_Go_enrichment

Down_Go_enrichment %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio,c("enrichedGene","inputGene"),"/") %>%
  tidyr::separate(BgRatio,c("pathBG","allBG"),"/") %>%
  dplyr::mutate(percent=100*as.numeric(enrichedGene)/as.numeric(pathBG)) %>%
  tidyr::separate(geneID,paste("gene",1:20,sep="_"),"/") %>%
  dplyr::select(-ONTOLOGY,-pvalue,-qvalue,-Count,-ID,-enrichedGene,-inputGene,-pathBG,-allBG) %>%
  tidyr::gather(-Description,-percent,-p.adjust,key="title",value="gene_id.x") %>%
  tidyr::drop_na() %>%
  dplyr::select(-title) %>%
  dplyr::inner_join(EHZ2_CBX2_common_targets.DE,by="gene_id.x") %>%
  dplyr::mutate(Description_p=paste(Description," (",round(percent,2),"%)",sep="")) -> Down_Go_enrichment.plotready
Down_Go_enrichment %>%
  as.data.frame() %>%
  dplyr::arrange(p.adjust) %>%
  head(20) %>% 
  dplyr::left_join(Down_Go_enrichment.plotready,by="Description") %>%
  dplyr::select(Description,Description_p) %>%
  unique() -> Down_Go_enrichment.up20
Down_Go_enrichment.plotready %>%
  dplyr::filter(Description %in% Down_Go_enrichment.up20$Description) %>%
  dplyr::group_by(gene_id.x) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::select(gene_id.x,n) %>% unique() -> generank
generank %>%
  dplyr::filter(gene_id.x %in% EZH2_CBX2_common_peak_common_targets)
Down_Go_enrichment.plotready %>%
  dplyr::filter(Description %in% Down_Go_enrichment.up20$Description) %>%
  ggplot(aes(x=gene_id.x,y=Description_p)) +
  geom_tile(aes(fill = log2FC),colour = "white") +
  scale_y_discrete(limits=Down_Go_enrichment.up20$Description_p) +
  scale_x_discrete(limits=generank$gene_id.x) +
  scale_fill_gradientn(colours=c("#00BFFF","red"),
                       name = "log2FC",
                       breaks=c(-2,0,2,4)) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        panel.background = element_blank(),
        # panel.grid.minor =element_line(color = "black",size=0.2),
        # panel.grid = element_line(colour = "grey", linetype = "dashed"),
        # panel.grid.major = element_line(
        #   colour = "grey",
        #   linetype = "dashed",
        #   size = 0.2
        # ),
        panel.border =element_rect(fill='transparent', color='black'))
ggsave(file.path(out_path,"Figure5","DOWN_EZH2-CBX2-targetgene-GOenrichment.heatmap.Rplot.pdf"),width = 12,height = 4)

##  KEGG enrichment ----
EHZ2_CBX2_common_targets.DE_info %>%
  tidyr::drop_na() %>%
  dplyr::filter(prob>0.9) -> EHZ2_CBX2_common_targets.prob0.9
EHZ2_CBX2_common_targets.prob0.9$log2FC ->x
names(x) = EHZ2_CBX2_common_targets.prob0.9$entrez_id
sort(x,decreasing = TRUE) ->x

kk <- gseKEGG(x, nPerm=1000,pvalueCutoff=1)
kk %>% 
  as.data.frame() %>%
  tidyr::separate(core_enrichment,paste("gene",1:100,sep="_"),"/") %>%
  dplyr::select(-ID,-setSize,-NES,-pvalue,-rank,-leading_edge) %>%
  tidyr::gather(-Description,-enrichmentScore,-p.adjust,-qvalues,key="title",value="entrez_id") %>%
  tidyr::drop_na() %>%
  dplyr::select(-title) %>%
  dplyr::inner_join(EHZ2_CBX2_common_targets.prob0.9,by="entrez_id") -> kk_info
kk_info %>%
  dplyr::filter(p.adjust<=0.05) %>% # enrichment showed p significant
  dplyr::select(Description,enrichmentScore,p.adjust,gene_id.x,log2FC) %>%
  unique() %>%
  dplyr::mutate(color=ifelse(log2FC>0,"red","blue"))-> kk_plotready
