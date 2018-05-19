data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets"
data_path_2 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"

# laod pac ----------------------------------------------------------------

library(magrittr)
library(clusterProfiler)
# load data ---------------------------------------------------------------

TF <- readr::read_tsv(file.path(data_path,"common-targets-180426-new","DOWN_TF_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::select(gene_id.x,entrez_id,log2FC) 
progene <- readr::read_tsv(file.path(data_path,"common-targets-180426-new","DOWN_pro_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::select(gene_id.x,entrez_id,log2FC) 
TF %>%
  rbind(progene) -> all_gene_de_info

# class genes into up and down part ---------------------------------------

library(org.Hs.eg.db)
library(clusterProfiler)
# all_gene.id <- bitr(c(all_gene_de_info$gene_id), fromType = "ALIAS",
#                    toType = c("ENTREZID"),
#                    OrgDb = org.Hs.eg.db)
# all_gene_de_info %>%
#   dplyr::rename("ALIAS"="gene_id") %>%
#   dplyr::inner_join(all_gene.id,by="ALIAS") %>%
#   as.data.frame()-> all_gene_info


# rownames(up_gene_info) <- up_gene_info$ENTREZID
# rownames(down_gene_info) <- down_gene_info$ENTREZID
# up_gene_info[,2] %>%
#   t() -> up_gene_info_filter
# colnames(up_gene_info_filter) <- up_gene_info$ENTREZID
# down_gene_info[,2] %>%
#   t() -> down_gene_info_filter
# colnames(down_gene_info_filter) <- down_gene_info$ENTREZID

# enrichment --------------------------------------------------------------

# GO -------------------------------------------------------------------
ego <- enrichGO(all_gene_de_info$entrez_id, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
library(enrichplot)
# goplot(ego)
# barplot(ego, showCategory=20)
dotplot(ego, showCategory=30)

go <- enrichGO(all_gene_de_info$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
library(ggplot2)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")


## remove redundent GO terms
# ego2 <- simplify(go)
# cnetplot(ego2, foldChange=geneList)
# cnetplot(ego2, foldChange=NULL, circular = TRUE, colorEdge = TRUE)
# upsetplot(ego)

heatplot(go, foldChange=all_gene_info$log2FC)
# emapplot(ego2)
# david <- enrichDAVID(gene = up_gene_info$ENTREZID,
#                      idType = "ENTREZ_GENE_ID",
#                      annotation = "KEGG_PATHWAY",
#                      david.user = "clusterProfiler@hku.hk")

