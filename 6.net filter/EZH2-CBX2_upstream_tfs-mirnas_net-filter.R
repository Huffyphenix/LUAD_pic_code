FFL_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/FFL/EZH2_CBX2_upstream"
# TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")

# load data ---------------------------------------------------------------
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute <- readr::read_tsv(file.path(FFL_data_path,"attribute.txt"),col_names = F)
network <- readr::read_tsv(file.path(FFL_data_path,"network.txt"),col_names = F) 

network %>% 
  dplyr::rename("source"="X1","target"="X2","regulate_type"="X3") %>%
  readr::write_tsv(file.path(FFL_data_path,"network.txt"))

TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
TSG <- readr::read_tsv(file.path(TSG_onco_data_path,"TSG.source_clear(at least two evidence-no confuse).tsv")) %>%
  dplyr::mutate(hallmark="TSG")
oncogene <- readr::read_tsv(file.path(TSG_onco_data_path,"oncogene.source_clear(at least two evidence-no confuse).tsv")) %>%
  dplyr::mutate(hallmark="oncogene")
confuse_gene <- readr::read_tsv(file.path(TSG_onco_data_path,"confuse_gene.source_clear(at least two evidence).tsv")) %>%
  dplyr::mutate(hallmark="confuse") %>%
  dplyr::select(symbol,hallmark) %>%
  dplyr::rename("SYMBOL"="symbol")
TSG %>%
  rbind(oncogene) %>%
  dplyr::select(SYMBOL,hallmark) %>%
  rbind(confuse_gene) -> all_cancer_relate_genes_noconfused

# TSGene <- readr::read_tsv(file.path(TSG_onco_data_path,"TSGene","LUAD_downregulate_TSG_in_TCGA-TSGeneDatabase.txt"))  %>%
#   dplyr::mutate(source="TSGene") %>%
#   dplyr::rename("symbol"="GeneSymbol","geneID"="GeneID") %>%
#   dplyr::select(symbol,geneID,source) 
# oncogene_database <- readr::read_tsv(file.path(TSG_onco_data_path,"oncogene database","oncogene_database-all_the_human_oncogenes.txt")) %>%
#   dplyr::mutate(source="oncogene_database") %>%
#   dplyr::rename("symbol"="OncogeneName","geneID"="OncogeneID") %>%
#   dplyr::select(symbol,geneID,source) 
# TSGene %>%
#   dplyr::inner_join(oncogene_database,by="geneID") %>%
#   .$geneID -> confused_gene
# TSGene %>%
#   dplyr::filter(! geneID %in% confused_gene) -> TSGene_noconfuse
# oncogene_database %>%
#   dplyr::filter(! geneID %in% confused_gene) -> oncogene_noconfuse
# TSGene_noconfuse %>%
#   rbind(oncogene_noconfuse) -> all_cancer_relate_genes_noconfused

# add hallmark info -------------------------------------------------------

attribute %>%
  dplyr::rename("SYMBOL"="X1","gene_type"="X2") %>%
  dplyr::left_join(all_cancer_relate_genes_noconfused,by="SYMBOL") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::mutate(mark=ifelse(hallmark=="TSGene",2,1)) %>%
  dplyr::mutate(mark=ifelse(is.na(hallmark),0,mark)) %>%
  dplyr::select(SYMBOL,gene_type,log2FC,mark) -> attribute.hallmark
attribute.hallmark %>%
  readr::write_tsv(file.path(FFL_data_path,"attribute.multi-source.hallmark-added.txt"))
