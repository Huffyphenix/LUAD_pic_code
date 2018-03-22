# data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-dwon_pro-TFgene"
data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-dwon_pro-TFgene"
data_path_2 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
data_path_3 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
# laod pac ----------------------------------------------------------------

library(magrittr)
library(clusterProfiler)
# load data ---------------------------------------------------------------
# DE gene----
TF <- readr::read_tsv(file.path(data_path,"NOIseq_DE_TF_FC2_cpm30.txt"))
progene <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_ProGene_FC2_cpm_30")) %>%
  dplyr::select(Gene_id,log2FC) %>%
  dplyr::rename("gene_id"="Gene_id")
TF %>%
  rbind(progene) -> all_gene_de_info

# all gene no filter: cpm>1
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil
# class genes into up and down part ---------------------------------------
all_gene_nofil %>%
  dplyr::filter(log2FC!="NA") %>%
  dplyr::filter(prob>0.9) -> all_gene_prob0.9

all_gene_de_info %>% 
  dplyr::filter(log2FC>0) -> up_gene_info

all_gene_de_info %>%
  dplyr::filter(log2FC<0) -> down_gene_info

library(org.Hs.eg.db)
up_gene.id <- bitr(c(up_gene_info$gene_id), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
down_gene.id <- bitr(c(down_gene_info$gene_id), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)
all_gene_prob0.9.id <- bitr(c(all_gene_prob0.9$gene_id), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)

up_gene_info %>%
  dplyr::rename("SYMBOL"="gene_id") %>%
  dplyr::inner_join(up_gene.id,by="SYMBOL") %>%
  as.data.frame()-> up_gene_info

down_gene_info %>%
  dplyr::rename("SYMBOL"="gene_id") %>%
  dplyr::inner_join(down_gene.id,by="SYMBOL") %>%
  as.data.frame() -> down_gene_info

all_gene_prob0.9 %>%
  dplyr::rename("SYMBOL"="gene_id") %>%
  dplyr::inner_join(all_gene_prob0.9.id,by="SYMBOL") %>%
  as.data.frame() -> all_gene_prob0.9_info
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
ego <- enrichGO(up_gene_info$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
library(enrichplot)
# goplot(ego)
# barplot(ego, showCategory=20)
dotplot(ego, showCategory=30)

go <- enrichGO(up_gene_info$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
library(ggplot2)
dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## down gene ---
ego_down <- enrichGO(down_gene_info$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
ego_down_si <- simplify(ego_down)
dotplot(ego_down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## remove redundent GO terms
# ego2 <- simplify(ego)
# cnetplot(ego2, foldChange=geneList)
# cnetplot(ego2, foldChange=NULL, circular = TRUE, colorEdge = TRUE)
# upsetplot(ego)
# heatplot(ego2, foldChange=NULL)
# emapplot(ego2)
# david <- enrichDAVID(gene = up_gene_info$ENTREZID,
#                      idType = "ENTREZ_GENE_ID",
#                      annotation = "KEGG_PATHWAY",
#                      david.user = "clusterProfiler@hku.hk")
all_gene_prob0.9_info %>%
  dplyr::filter(abs(log2FC)>1) -> all_gene_prob0.9fc2_info

all_gene_prob0.9fc2_info[,6] ->x
names(x) = all_gene_prob0.9fc2_info$ENTREZID
sort(x,decreasing = TRUE) ->x

kk <- gseKEGG(x, nPerm=1000)
ridgeplot(kk,showCategory = 30)
gseaplot(kk, geneSetID = 1, title = kk$Description[1])
heatplot(kk)
ego3 <- gseGO(geneList     = x,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 10000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ridgeplot(ego3)
