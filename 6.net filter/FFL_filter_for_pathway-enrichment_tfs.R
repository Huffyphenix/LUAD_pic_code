FFL_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL"
TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

# load data ---------------------------------------------------------------
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute <- read.table(file.path(FFL_data_path,"attribute-net-all.txt")) %>%
  dplyr::rename("SYMBOL"="V1","gene_type"="V2","up_down"="V3") %>%
  dplyr::as.tbl()
network <- read.table(file.path(FFL_data_path,"network.txt")) %>%
  dplyr::rename("From"="V1","To"="V2","regulate_type"="V3") %>%
  dplyr::as.tbl()
network %>% 
  readr::write_tsv(file.path(FFL_data_path,"network.txt"))
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
  rbind(confuse_gene) -> all_cancer_relate_genes


# filter ------------------------------------------------------------------
attribute %>%
  dplyr::left_join(all_cancer_relate_genes,by="SYMBOL") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,gene_type,hallmark,log2FC) -> attribute.hallmark

attribute.hallmark %>%
  readr::write_tsv(file.path(FFL_data_path,"attribute.hallmark-added.txt"))
                   