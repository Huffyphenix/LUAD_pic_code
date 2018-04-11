# data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-dwon_pro-TFgene"
data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-dwon_pro-TFgene"
data_path_1<- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
data_path_2 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
# laod pac ----------------------------------------------------------------

library(magrittr)
library(clusterProfiler)
# load data ---------------------------------------------------------------
# DE gene----
TF <- readr::read_tsv(file.path(data_path_2,"NOISeq_DE_TF_FC2_cpm_30")) %>%
  dplyr::select(gene_id,log2FC) 
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

all_gene_prob0.9.id <- bitr(c(all_gene_prob0.9$gene_id), fromType = "ALIAS",
                     toType = c("SYMBOL"),
                     OrgDb = org.Hs.eg.db) %>%
  .$SYMBOL %>% unique() %>%
  bitr(fromType = "SYMBOL",
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
  dplyr::filter(abs(log2FC)>=1) -> all_gene_prob0.9fc2_info

# only for fc 2 ----
all_gene_prob0.9fc2_info[,6] ->x
names(x) = all_gene_prob0.9fc2_info$ENTREZID
sort(x,decreasing = TRUE) ->x

kk_fc2 <- gseKEGG(x, nPerm=1000,pvalueCutoff=1)
kk_fc2 %>% 
  as.data.frame() %>%
  tidyr::separate(core_enrichment,paste("gene",1:100,sep="_"),"/") %>%
  dplyr::select(-ID,-setSize,-NES,-pvalue,-rank,-leading_edge) %>%
  tidyr::gather(-Description,-enrichmentScore,-p.adjust,-qvalues,key="title",value="ENTREZID") %>%
  tidyr::drop_na() %>%
  dplyr::select(-title) %>%
  dplyr::inner_join(all_gene_prob0.9_info,by="ENTREZID") -> kk_fc2_info
kk_fc2_info %>%
  dplyr::filter(p.adjust<=0.05) %>%
  dplyr::select(Description,enrichmentScore,p.adjust,SYMBOL,log2FC) %>%
  unique() %>%
  dplyr::mutate(color=ifelse(log2FC>0,"red","blue"))-> kk_fc2_plotready
data_path_4 <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-FC2-prob0.9-kegg-gsea"
kk_fc2_plotready$SYMBOL %>%
  bitr(fromType = "SYMBOL",
       toType = c("UNIPROT"),
       OrgDb = org.Hs.eg.db) %>%
  dplyr::inner_join(kk_nofc_plotready,by="SYMBOL") %>%
  readr::write_tsv(file.path(data_path_4,"kk_fc2_uniprot_color_for_keggmapper_padjust0.05.tsv"))


kk_fc2_plotready %>%
  dplyr::arrange(enrichmentScore) %>%
  dplyr::select(Description) %>%
  unique() -> kk_fc2_description_rank
readr::write_tsv(kk_nofc_description_rank,file.path(data_path_1,"kk_nofc_description_rank_padjust0.05"))
kk_nofc_rank <- readr::read_tsv(file.path(data_path_1,"kk_nofc_description_rank_padjust0.05"))

library(ggplot2)
kk_fc2_plotready %>%
  ggplot(aes(x=log2FC,y=Description)) +
  ggridges::geom_density_ridges_gradient(aes(fill = enrichmentScore), scale = 3, size = 0.1,rel_min_height = 0.01) +
  scale_fill_gradientn(colours=c("#00BFFF","white","red"),name = "enrichmentScore") +
  scale_y_discrete(limits=kk_fc2_description_rank$Description) +
  ylab("KEGG pathway") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        ),
        panel.border =element_rect(fill='transparent', color='black'))
ggsave(file.path(data_path_4,"gseaKEGG_for_LUAD-FC2-prob0.9_padjust0.05.pdf"),width = 8,height = 8)

kk_fc2_plotready %>%
  dplyr::filter(enrichmentScore>0) %>%
  ggplot(aes(x=SYMBOL,y=Description)) +
  geom_tile(aes(fill = log2FC)) +
  scale_y_discrete(limits=c("Cell cycle",
                            "p53 signaling pathway","Oocyte meiosis",
                            "Progesterone-mediated oocyte maturation",
                            "Homologous recombination",
                            "Fanconi anemia pathway")) +
  scale_fill_gradientn(colours=c("#00BFFF","red"),
                       name = "log2FC",
                       breaks=c(0,2,4,6,8)) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        ),
        panel.border =element_rect(fill='transparent', color='black'))
ggsave(file.path(data_path_4,"heatmap_gseaKEGG_for_LUAD-FC2-prob0.9_padjust0.05_upregulate.pdf"),width = 10,height = 3)
# ridgeplot(kk_fc2,showCategory = 30)
# gseaplot(kk, geneSetID = 3, title = kk$Description[3])
# heatplot(kk,foldChange = x)

# for all prob 0.9 ----
all_gene_prob0.9_info[,6] ->y
names(y) = all_gene_prob0.9_info$ENTREZID
sort(y,decreasing = TRUE) ->y

kk_nofc <- gseKEGG(y, nPerm=1000,pvalueCutoff=1)
ridgeplot(kk_nofc,showCategory = 40)
kk_nofc %>% 
  as.data.frame() %>%
  tidyr::separate(core_enrichment,paste("gene",1:100,sep="_"),"/") %>%
  dplyr::select(-ID,-setSize,-NES,-pvalue,-rank,-leading_edge) %>%
  tidyr::gather(-Description,-enrichmentScore,-p.adjust,-qvalues,key="title",value="ENTREZID") %>%
  tidyr::drop_na() %>%
  dplyr::select(-title) %>%
  dplyr::inner_join(all_gene_prob0.9_info,by="ENTREZID") -> kk_nofc_info
library(RColorBrewer)
Colormap<- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
kk_nofc_info %>%
  dplyr::filter(p.adjust<=0.05) %>%
  dplyr::select(Description,enrichmentScore,p.adjust,SYMBOL,log2FC) %>%
  unique() %>%
  dplyr::mutate(color=ifelse(log2FC>0,"red","blue"))-> kk_nofc_plotready
kk_nofc_plotready$SYMBOL %>%
  bitr(fromType = "SYMBOL",
       toType = c("UNIPROT"),
       OrgDb = org.Hs.eg.db) %>%
  dplyr::inner_join(kk_nofc_plotready,by="SYMBOL") %>%
  readr::write_tsv(file.path(data_path_1,"kk_nofc_uniprot_color_for_keggmapper_padjust0.05.tsv"))


kk_nofc_plotready %>%
  dplyr::arrange(enrichmentScore) %>%
  dplyr::select(Description) %>%
  unique() -> kk_nofc_description_rank
readr::write_tsv(kk_nofc_description_rank,file.path(data_path_1,"kk_nofc_description_rank_padjust0.05"))
kk_nofc_rank <- readr::read_tsv(file.path(data_path_1,"kk_nofc_description_rank_padjust0.05"))

library(ggplot2)
kk_nofc_plotready %>%
  ggplot(aes(x=log2FC,y=Description)) +
  ggridges::geom_density_ridges_gradient(aes(fill = enrichmentScore), scale = 3, size = 0.1,rel_min_height = 0.01) +
  scale_fill_gradientn(colours=c("#00BFFF","white","red"),name = "enrichmentScore") +
  scale_y_discrete(limits=kk_nofc_rank$Description) +
  ylab("KEGG pathway") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        ),
        panel.border =element_rect(fill='transparent', color='black'))
ggsave(file.path(data_path_1,"gseaKEGG_for_LUAD-noFC-prob0.9_padjust0.05.pdf"),width = 8,height = 8)

# heatmap 
kk_nofc_plotready %>%
  dplyr::filter(enrichmentScore>0) %>%
  ggplot(aes(x=SYMBOL,y=Description)) +
  geom_tile(aes(fill = log2FC)) +
  scale_fill_gradientn(colours=c("#00BFFF","red"),
                       name = "log2FC",
                       breaks=c(0,2,4,6,8))

gseaplot(kk, geneSetID = 3, title = kk$Description[3])
heatplot(kk,foldChange = y)

# network
kk_nofc_plotready %>%
  dplyr::filter(SYMBOL %in% TF$gene_id & enrichmentScore>0) %>%
  dplyr::select(Description,SYMBOL) -> TF_crosstalk
TF_crosstalk %>%
  dplyr::select(Description) %>%
  dplyr::mutate(type=1) %>%
  dplyr::rename("node"="Description") %>%
  unique() -> TF_crosstalk.pathway.node
TF_crosstalk %>%
  dplyr::select(SYMBOL) %>%
  dplyr::mutate(type=2) %>%
  dplyr::rename("node"="SYMBOL") %>%
  unique() -> TF_crosstalk.TF.node
TF_crosstalk.pathway.node %>%
  rbind(TF_crosstalk.TF.node) -> TF_crosstalk.node
TF_crosstalk.node %>%
  readr::write_tsv(file.path(data_path_1,"TF.crosstalk","TF_crosstalk.node.attribute.txt"))
TF_crosstalk %>%
  readr::write_tsv(file.path(data_path_1,"TF.crosstalk","TF_crosstalk.node.network.txt"))

ego3 <- gseGO(geneList     = x,
              OrgDb        = org.Hs.eg.db,
              ont          = "all",
              nPerm        = 10000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ridgeplot(ego3)
