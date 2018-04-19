data_path_3 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
data_path <- "F:/?ҵļ?????/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
# load data ---------------------------------------------------------------
## DE genes
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil

TF <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_TF_FC2_cpm_30")) %>%
  dplyr::select(gene_id,log2FC) 
progene <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_ProGene_FC2_cpm_30")) %>%
  dplyr::select(Gene_id,log2FC) %>%
  dplyr::rename("gene_id"="Gene_id")
TF %>%
  rbind(progene) -> all_gene_de_info
all_gene_de_info %>% 
  dplyr::filter(log2FC>0) -> up_gene_info

all_gene_de_info %>%
  dplyr::filter(log2FC<0) -> down_gene_info
up_gene.id <- bitr(c(up_gene_info$gene_id), fromType = "ALIAS",
                   toType = c("SYMBOL"),
                   OrgDb = org.Hs.eg.db) %>%
  .$SYMBOL %>% unique() %>%
  bitr(fromType = "SYMBOL",
       toType = c("ENTREZID"),
       OrgDb = org.Hs.eg.db)
down_gene.id <- bitr(c(down_gene_info$gene_id), fromType = "ALIAS",
                                 toType = c("SYMBOL"),
                                 OrgDb = org.Hs.eg.db) %>%
  .$SYMBOL %>% unique() %>%
  bitr(fromType = "SYMBOL",
       toType = c("ENTREZID"),
       OrgDb = org.Hs.eg.db)
## onco and TSG info
TSGene <- readr::read_tsv(file.path(TSG_onco_data_path,"TSGene","all_tumor_supressor.txt"))  %>%
  dplyr::mutate(source="TSGene") %>%
  dplyr::rename("symbol"="GeneSymbol","geneID"="GeneID") %>%
  dplyr::select(symbol,geneID,source) 
TSGene_luad <- readr::read_tsv(file.path(TSG_onco_data_path,"TSGene","LUAD_downregulate_TSG_in_TCGA-TSGeneDatabase.txt"))  %>%
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
up_gene.id %>%
  dplyr::as.tbl() %>%
  dplyr::filter(ENTREZID %in% TSGene_noconfuse$geneID) %>%
  nrow() -> up_TSG
up_gene.id %>%
  dplyr::as.tbl() %>%
  dplyr::filter(ENTREZID %in% oncogene_noconfuse$geneID) %>%
  nrow() -> up_onco
down_gene.id %>%
  dplyr::as.tbl() %>%
  dplyr::filter(ENTREZID %in% TSGene_noconfuse$geneID) %>%
  nrow() -> down_TSG
down_gene.id %>%
  dplyr::as.tbl() %>%
  dplyr::filter(ENTREZID %in% oncogene_noconfuse$geneID) %>%
  nrow() -> down_onco

matrix(c(33,1333,10,1391),nrow = 2) -> x
chisq.test(x,correct = T)
