source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
library(pathview)

# data path ---------------------------------------------------------------

data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-dwon_pro-TFgene"
data_path_1<- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
data_path_3 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

# all gene no filter: cpm>1
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil

all_gene_nofil %>%
  dplyr::filter(log2FC!="NA") %>%
  dplyr::filter(prob>0.9) -> all_gene_prob0.9

library(org.Hs.eg.db)
library(clusterProfiler)
all_gene_prob0.9.id <- bitr(c(all_gene_prob0.9$gene_id), fromType = "ALIAS",
                            toType = c("SYMBOL"),
                            OrgDb = org.Hs.eg.db) %>%
  .$SYMBOL %>% unique() %>%
  bitr(fromType = "SYMBOL",
       toType = c("ENTREZID"),
       OrgDb = org.Hs.eg.db)
all_gene_prob0.9 %>%
  dplyr::rename("SYMBOL"="gene_id") %>%
  dplyr::inner_join(all_gene_prob0.9.id,by="SYMBOL") %>%
  as.data.frame() -> all_gene_prob0.9_info
all_gene_prob0.9_info %>%
  dplyr::select(log2FC) %>%
  as.matrix() -> data
rownames(data) <- all_gene_prob0.9_info$ENTREZID

# map04110 cell cycle -----------------------------------------------------
setwd("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure1")
pv.out <- pathview(gene.data = data, pathway.id = "04110",species = "hsa", out.suffix = "cellcycle_DEG",limit=list(gene=2))

# map04115 p53----------------------------------------------------------------
pv.out <- pathview(gene.data = data, pathway.id = "04115",species = "hsa", out.suffix = "p53_DEG",limit=list(gene=2))

