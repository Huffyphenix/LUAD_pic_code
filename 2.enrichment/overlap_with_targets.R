data_path_1 <- c("F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea")
data_path_2 <- c("F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets")
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
data_path_4 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

# load data ---------------------------------------------------------------

gseaKEGG_result_all <- readr::read_tsv(file.path(data_path_1,"gseaKEGG_result-gather.tsv"))
gseaKEGG_result <- readr::read_tsv(file.path(data_path_1,"gseaKEGG_results.tsv"))
EZH2_CBX2_common_Pro <- readr::read_tsv(file.path(data_path_2,"common-targets-180426-new","DOWN1.5_allexp30_pro_EHZ2_CBX2_common_targets.DE_info"))
EZH2_CBX2_common_TF <- readr::read_tsv(file.path(data_path_2,"common-targets-180426-new","DOWN1.5_allexp30_TF_EHZ2_CBX2_common_targets.DE_info"))
EZH2_CBX2_common_targets <- rbind(EZH2_CBX2_common_Pro,EZH2_CBX2_common_TF)
EZH2_targets <- readr::read_tsv(file.path(data_path_2,"common-targets-180426-new","EZH2_targets_proteincoding.DE_info"))
CBX2_targets <- readr::read_tsv(file.path(data_path_2,"common-targets-180426-new","CBX2_targets_proteincoding.DE_info"))

TF <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_FC2_cpm_30")) %>%
  dplyr::select(gene_id,log2FC) 
progene <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_FC2_cpm_30")) %>%
  dplyr::rename("gene_id"="Gene_id") %>%
  dplyr::select(gene_id,log2FC) 
TF %>%
  rbind(progene) -> all_gene_de_info

# all gene no filter: cpm>1
TF_nofil <- readr::read_tsv(file.path(data_path_4,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_4,"NOISeq_DE_ProGene_cpm_1_noFDR"))
rbind(TF_nofil,progene_nofil) -> all_gene_nofil

all_gene_nofil %>%
  dplyr::filter(log2FC!="NA") %>%
  dplyr::filter(prob>0.9) -> all_gene_prob0.9

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

# filter ------------------------------------------------------------------

interst_path <- "PPAR signaling pathway"

gseaKEGG_result_all %>%
  dplyr::filter(Description %in% interst_path) -> PPAR

PPAR %>%
  dplyr::filter(SYMBOL %in% EZH2_CBX2_common_targets$gene_id.y)

##### common targets in downregulate pathway --------
gseaKEGG_result_all %>%
  dplyr::filter(SYMBOL %in% EZH2_CBX2_common_targets$gene_id.y) %>%
  dplyr::filter(p.adjust<0.05 & enrichmentScore<0) %>%
  dplyr::arrange(Description) -> EZH2_CBX2_common_targets_in_down_kegg
EZH2_CBX2_common_targets_in_down_kegg %>%
  readr::write_tsv(file.path(data_path_2,"common-targets-180426-new/overlap gseaKEGG","EZH2_CBX2_common_targets_in_down_kegg.tsv"))

##### EZH2 targets in downregulate pathway --------
PPAR %>%
  dplyr::filter(SYMBOL %in% EZH2_targets$gene_id.y)

gseaKEGG_result_all %>%
  dplyr::filter(SYMBOL %in% EZH2_targets$gene_id.y) %>%
  dplyr::filter(p.adjust<0.05 & enrichmentScore<0) %>%
  dplyr::arrange(Description) -> EZH2_targets_in_down_kegg
EZH2_targets_in_down_kegg %>%
  readr::write_tsv(file.path(data_path_2,"common-targets-180426-new/overlap gseaKEGG","EZH2_targets_in_down_kegg.tsv"))

##### CBX2 targets in downregulate pathway --------
gseaKEGG_result_all %>%
  dplyr::filter(SYMBOL %in% CBX2_targets$gene_id.y) %>%
  dplyr::filter(p.adjust<0.05 & enrichmentScore<0) %>%
  dplyr::arrange(Description) -> CBX2_targets_in_down_kegg
CBX2_targets_in_down_kegg %>%
  readr::write_tsv(file.path(data_path_2,"common-targets-180426-new/overlap gseaKEGG","CBX2_targets_in_down_kegg.tsv"))

# pathview ----------------------------------------------------------------

.libPaths("H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(pathview)
setwd("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure1")
pv.out<- pathview(gene.data = data, pathway.id = "03320",species = "hsa", out.suffix = "PPAR_DEG",limit=list(gene=2))
