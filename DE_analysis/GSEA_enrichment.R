
# enrichment analysis -----------------------------------------------------
library(magrittr)
library(clusterProfiler)
# for gene DE between CBX2-EZH2 high and low

# path --------------------------------------------------------------------
# HUST
basic_path <- file.path("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD")
# HOME
basic_path <- file.path("F:/胡斐斐/我的坚果云/ENCODE-TCGA-LUAD/")

# exp_path <- "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
res_path <- file.path(basic_path,"Figure/DE_between_high_low_CBX2_EZH2")

# load data ---------------------------------------------------------------

CBX2_EZH2_DE <- readr::read_tsv(file.path(res_path,"CBX2_and_EZH2","DE_gene_res.tsv"))

# GSEA enrichment ---------------------------------------------------------
# p < 0.05
CBX2_EZH2_DE %>%
  dplyr::filter(p.value<=0.05) -> CBX2_EZH2_DE.p.0.05

y <- CBX2_EZH2_DE.p.0.05$log2FC_HvsL
names(y) <- CBX2_EZH2_DE.p.0.05$entrez_id
sort(y,decreasing = TRUE) ->y # 20530 genes

kk_nofc <- gseKEGG(y, nPerm=1000,pvalueCutoff=1)  # pvaluecutoff is important for result reproducible. must 0.05
ridgeplot(kk_nofc,showCategory = 40)
kk_nofc %>% 
  as.data.frame() %>%
  tidyr::separate(core_enrichment,paste("gene",1:100,sep="_"),"/") %>%
  dplyr::select(-ID,-setSize,-NES,-pvalue,-rank,-leading_edge) %>%
  tidyr::gather(-Description,-enrichmentScore,-p.adjust,-qvalues,key="title",value="entrez_id") %>%
  tidyr::drop_na() %>%
  # dplyr::select(-title) %>%
  dplyr::inner_join(CBX2_EZH2_DE %>% dplyr::mutate(entrez_id=as.character(entrez_id)),by="entrez_id") -> kk_nofc_info
kk_nofc %>% 
  as.data.frame() -> kk_nofc.tible
kk_nofc_info %>%
  dplyr::select(Description,title,symbol) %>%
  tidyr::spread(key = title,value=symbol) %>%
  tidyr::unite("enriched_genes",paste("gene",1:86,sep="_"),sep="/") %>%
  dplyr::mutate(enriched_genes=gsub("/NA","",enriched_genes)) %>%
  dplyr::inner_join(kk_nofc.tible,by="Description") %>%
  dplyr::mutate(Counts = strsplit(enriched_genes,"/") %>%  lapply(length) %>% unlist()) %>%
  dplyr::mutate(`Percent (%)`=round(Counts*100/setSize,0)) %>%
  dplyr::select(Description,ID,setSize,enrichmentScore,enriched_genes,Counts,`Percent (%)`,pvalue,p.adjust,qvalues) -> kk_nofc_info.1

kk_nofc_info.1 %>%
  dplyr::filter(pvalue<=0.05) %>%
  dplyr::filter(enrichmentScore<0) %>%
  unique() %>%
  dplyr::arrange(enrichmentScore) -> kk_nofc_plotready.down

gseaplot(kk_nofc,geneSetID = "hsa03320",by = "all", title = "PPAR signaling pathway in tumor with\nCBX2 and EZH2 both higher; p = 0.02, q = 0.1")

gseaplot(kk_nofc,geneSetID = "hsa04110",by = "all", title = "Cell cycle in tumor with CBX2 and EZH2 both higher;\np = 0.001, q = 0.03")

library(org.Hs.eg.db)
mapped <- mappedkeys(org.Hs.egPATH2EG)
L <- as.list(org.Hs.egPATH2EG[mapped])
L[["03320"]] -> pparg_genes

CBX2_EZH2_DE.p.0.05 %>%
  dplyr::arrange(desc(log2FC_HvsL)) %>%
  dplyr::mutate(rank = 1:12226,entrez_id=as.character(entrez_id)) %>%
  # dplyr::select(entrez_id,rank) %>%
  dplyr::filter(entrez_id %in% pparg_genes) -> pparg_genes_rank


kk_nofc_info %>%
  dplyr::filter(Description == "PPAR signaling pathway") %>%
  dplyr::inner_join(pparg_genes_rank, by = "entrez_id")

# <0.05 & downregulate
CBX2_EZH2_DE %>%
  dplyr::filter(p.value<=0.05 & log2FC_HvsL <0) -> CBX2_EZH2_DE.p.0.05.down

y <- CBX2_EZH2_DE.p.0.05.down$log2FC_HvsL
names(y) <- CBX2_EZH2_DE.p.0.05.down$entrez_id
sort(y,decreasing = TRUE) ->y # 20530 genes

kk_nofc <- gseKEGG(y, nPerm=1000,pvalueCutoff=1)  # pvaluecutoff is important for result reproducible. must 0.05
gseaplot(kk_nofc,geneSetID = "hsa03320",by = "all", title = "PPAR signaling pathway")
ridgeplot(kk_nofc,showCategory = 40)
