# data path ---------------------------------------------------------------

genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
data_path<- "H:/data"
library(magrittr)
# data manage -------------------------------------------------------------
genelist <- readr::read_tsv(file.path(genelist_path,"all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(prob>=0.99 & abs(log2FC)>=0.585)

Animal_TF <-  readr::read_tsv(file.path(data_path,"AnimalTFDB","Homo_sapiens_transcription_factors_gene_list.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls")) %>%
  .$Symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Animal_TF %>%
  dplyr::filter(! Entrez_ID %in% enzyme_lsit$ENTREZID) -> Animal_TF
Animal_TF$Entrez_ID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db) -> Animal_TF_entrez

genelist %>%
  dplyr::filter(entrez_id %in% Animal_TF_entrez$ENTREZID) %>%
  dplyr::filter(log2FC<0) -> DOWN_commone_targets.TF 
DOWN_commone_targets.TF %>%
  readr::write_tsv(file.path(genelist_path,"FFL/input","DOWN_commone_targets.TF"))
genelist %>%
  dplyr::filter(! entrez_id %in% Animal_TF_entrez$ENTREZID) %>%
  dplyr::filter(log2FC<0) -> DOWN_commone_targets.pro
DOWN_commone_targets.pro %>%
  readr::write_tsv(file.path(genelist_path,"FFL/input","DOWN_commone_targets.pro"))

# TSG and oncogene statistic ----------------------------------------------
TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
TSGene <- readr::read_tsv(file.path(TSG_onco_data_path,"TSGene","LUAD_downregulate_TSG_in_TCGA-TSGeneDatabase.txt"))  %>%
  dplyr::mutate(source="TSGene") %>%
  dplyr::rename("symbol"="GeneSymbol","geneID"="GeneID") %>%
  dplyr::select(symbol,geneID,source) 
oncogene_database <- readr::read_tsv(file.path(TSG_onco_data_path,"oncogene database","oncogene_database-all_the_human_oncogenes.txt")) %>%
  dplyr::mutate(source="oncogene_database") %>%
  dplyr::rename("symbol"="OncogeneName","geneID"="OncogeneID") %>%
  dplyr::select(symbol,geneID,source) 

TSGene %>%
  dplyr::inner_join(oncogene_database,by="geneID") %>%
  .$geneID -> confused_gene
TSGene %>%
  dplyr::filter(! geneID %in% confused_gene) -> TSGene_noconfuse
oncogene_database %>%
  dplyr::filter(! geneID %in% confused_gene) -> oncogene_noconfuse
TSGene_noconfuse %>%
  rbind(oncogene_noconfuse) -> all_cancer_relate_genes_noconfused

genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(gene_id.x %in% TSGene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(gene_id.x %in% oncogene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC>0) %>%
  dplyr::filter(gene_id.x %in% oncogene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC>0) %>%
  dplyr::filter(gene_id.x %in% TSGene_noconfuse$symbol) %>% nrow()


### >>>> construc FFL in server 1:/home/huff/LUAD_cancer/FFL_quantification_data/EZH2_CBX2_targets/
ffl_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets"
TF2miRNA <- readr::read_tsv(file.path(ffl_path,"TF2miRNA"),col_names = F)
TF2gene <- readr::read_tsv(file.path(ffl_path,"TF2gene"),col_names = F)
miRNA2TF <- readr::read_tsv(file.path(ffl_path,"miRNA2TF"),col_names = F)
miRNA2gene <- readr::read_tsv(file.path(ffl_path,"miRNA2gene"),col_names = F)
attribute <- readr::read_tsv(file.path(ffl_path,"attribute.txt"),col_names = F)

# filter ------------------------------------------------------------------

# only remain E2F1 and SOX4 as upstream of miRNAs, cause other TFs are all downregulate because of EZH2 and CBX2.
TF2miRNA %>%
  dplyr::filter(X1 %in% c("E2F1","SOX4")) %>%
  dplyr::mutate(X4=1) -> E2F1_SOX4_2_miRNA

# only remain miRNAs which regulate by E2F1 and SOX4

miRNA2TF %>%
  dplyr::filter(X1 %in% unique(E2F1_SOX4_2_miRNA$X2)) %>%
  dplyr::mutate(X3="notshow") %>%
  dplyr::mutate(X4=3) -> miRNA2TF.filter

miRNA2gene %>%
  dplyr::filter(X1 %in% unique(E2F1_SOX4_2_miRNA$X2)) %>%
  dplyr::mutate(X3="notshow") %>%
  dplyr::mutate(X4=4) -> miRNA2gene.filter

# only remain miRNAs which regulate by E2F1 and SOX4
TF2gene %>%
  dplyr::filter(! X1 %in% c("E2F1","SOX4")) %>%
  dplyr::filter(X2 %in% unique(miRNA2gene.filter$X2)) %>%
  dplyr::mutate(X4=2)-> TF2gene.filter

# combine -----------------------------------------------------------------

E2F1_SOX4_2_miRNA %>%
  rbind(miRNA2TF.filter) %>%
  rbind(miRNA2gene.filter) %>%
  rbind(TF2gene.filter) -> network

data.frame(gene=E2F1_SOX4_2_miRNA$X1,type=1) %>%
  rbind(data.frame(gene=E2F1_SOX4_2_miRNA$X2,type=2)) %>%
  rbind(data.frame(gene=miRNA2TF.filter$X2,type=1)) %>%
  rbind(data.frame(gene=miRNA2TF.filter$X1,type=2)) %>%
  rbind(data.frame(gene=miRNA2gene.filter$X1,type=2)) %>%
  rbind(data.frame(gene=miRNA2gene.filter$X2,type=3)) %>%
  rbind(data.frame(gene=TF2gene.filter$X1,type=1)) %>%
  rbind(data.frame(gene=TF2gene.filter$X2,type=3)) %>%
  dplyr::mutate(gene=as.character(gene)) %>%
  unique()  -> attribute


# add fold change info ----------------------------------------------------
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute %>%
  dplyr::rename("SYMBOL"="gene") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,type,log2FC) -> attribute.fc


# add TSG info ------------------------------------------------------------



attribute.fc %>%
  dplyr::rename("symbol"="SYMBOL") %>%
  dplyr::left_join(all_cancer_relate_genes_noconfused,by="symbol") %>%
  dplyr::mutate(mark=ifelse(source=="TSGene",2,1)) %>%
  dplyr::mutate(mark=ifelse(is.na(source),0,mark)) %>%
  dplyr::select(symbol,type,log2FC,mark) -> attribute.fc.tsg
  

attribute.fc.tsg %>%
  readr::write_tsv(file.path(genelist_path,"FFL/EZH2_CBX2_targets/√network.filter","attribute_fc_tsg.txt"))
network %>%
  readr::write_tsv(file.path(genelist_path,"FFL/EZH2_CBX2_targets/√network.filter","network.txt"))

# genes out of control of CNV, methylation, and miRNA regulation -----
# 1 methylation
mthy_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"
methy_regu_gene <- readr::read_tsv(file.path(mthy_path,"Figure4/Figure5","genes_regulate_by_methy.tsv"))

# 2 CNV 
Down_common_targets.CNV <- readr::read_tsv(file.path(data_path,"Figure4/Figure5","Down_common_targets.CNV-percent.tsv"))
Down_common_targets.CNV %>%
  dplyr::filter(homo_del>5) -> Down_common_targets.CNV_dele_5

# 3 miRNA regulation
attribute %>%
  dplyr::filter(type %in% c(1,3)) %>%
  dplyr::filter(! gene %in% c("E2F1","SOX4")) -> miRNA_regu_genes

# 4 filter
genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(! gene_id.x %in% methy_regu_gene$symbol) %>%
  dplyr::filter(! gene_id.x %in% Down_common_targets.CNV_dele_5$SAMPLE_ID) %>%
  dplyr::filter(! gene_id.x %in% miRNA_regu_genes$gene) %>%
  dplyr::filter(gene_id.x %in% TSGene_noconfuse$symbol)
