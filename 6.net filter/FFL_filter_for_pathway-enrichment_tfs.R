# E zHou ----
FFL_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-ppar-relatedgenes"
FFL_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-cellcycle-relatedgenes"

TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
enrich_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

#  HUST ----
FFL_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-ppar-relatedgenes"
FFL_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-cellcycle-relatedgenes"


TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_3 <- "G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
enrich_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
miRNA_exp_path <- "G:/data/TCGA/TCGA_data"

# FFL_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-ppar-relatedgenes"
# TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
# data_path_3 <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
# enrich_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
# .libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")

# load data ---------------------------------------------------------------
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute <- read.table(file.path(FFL_data_path,"attribute.txt")) %>%
  dplyr::rename("SYMBOL"="V1","gene_type"="V2") %>%
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



# load expression 
miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 


TSG %>%
  rbind(oncogene) %>%
  dplyr::select(SYMBOL,hallmark) %>%
  rbind(confuse_gene) -> all_cancer_relate_genes

cell_cycle_relate <- c("Cell cycle","Oocyte meiosis","DNA replication",
                       "Homologous recombination","p53 signaling pathway",
                       "Progesterone-mediated oocyte maturation")
ppar_relate <- c("PPAR signaling pathway")
enrichment<- readr::read_tsv(file.path(enrich_path,"gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description %in% ppar_relate)

# filter ------------------------------------------------------------------
attribute %>%
  dplyr::left_join(all_cancer_relate_genes,by="SYMBOL") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,gene_type,hallmark,log2FC) %>%
  dplyr::mutate(mark=ifelse(hallmark=="TSG",2,1)) %>%
  dplyr::mutate(mark=ifelse(is.na(hallmark),0,mark)) -> attribute.hallmark
library(plyr)
attribute %>%
  dplyr::left_join(enrichment,by="SYMBOL") %>%
  dplyr::select(SYMBOL,Description) %>%
  dplyr::mutate(Description=ifelse(is.na(Description),"NA",Description)) %>%
  unique() %>%
  dplyr::arrange(Description) %>%
  ddply(.(SYMBOL), summarise,
        Description=paste(Description,collapse=",")) %>%
  readr::write_tsv(file.path(FFL_data_path,"attribute.enrichment.txt"))
  
attribute.hallmark %>%
  readr::write_tsv(file.path(FFL_data_path,"attribute.hallmark-added.txt"))



# only for mirna regulation -----------------------------------------------
mirna_ppar <- readr::read_tsv(file.path(FFL_data_path,"miRNA2gene"),col_names = F)

mirna_ppar %>%
  tidyr::gather(key="type",value="gene") %>%
  dplyr::mutate(type=ifelse(type=="X1",1,2)) %>%
  unique() -> mirna_ppar_attribute
mirna_ppar %>%
  dplyr::rename("regulator" = "X1", "gene" = "X2") -> mirna_ppar_network

mirna_ppar_attribute %>%
  dplyr::rename("SYMBOL" = "gene") %>%
  dplyr::left_join(all_cancer_relate_genes,by="SYMBOL") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::rename("gene" = "SYMBOL") %>%
  dplyr::select(type,gene,hallmark,log2FC) %>%
  readr::write_tsv(file.path(FFL_data_path,"miRNA_PPAR_network","mirna_ppar_attribute.txt"))
mirna_ppar_network %>%
  dplyr::rename("SYMBOL" = "gene") %>%
  dplyr::left_join(all_cancer_relate_genes,by="SYMBOL") %>%
  dplyr::rename("gene" = "SYMBOL") %>%
  readr::write_tsv(file.path(FFL_data_path,"miRNA_PPAR_network","mirna_ppar_network.txt"))

# get expression 
miRNA_exp %>%
  dplyr::filter(name %in% c(mirna_ppar$X1)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> net_miRNA_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(mirna_ppar$X2)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_gene_exp.gather

# function for mirna2gene
fn_get_spm_a <- function(m_data,g_exp){
  g_exp %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(spm=purrr::map(data,m_data=m_data,fn_get_spm_b)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::ungroup() %>%
    dplyr::rename(Cor=estimate)-> .out
  return(.out)
}

fn_get_spm_b <- function(g_data,m_data){
  g_data %>%
    dplyr::inner_join(m_data,by="sample") ->tmp
  tmp %>%
    dplyr::mutate(gene_exp=gene_exp+runif(nrow(tmp),min=0,max=0.001)) %>%
    dplyr::mutate(mirna_exp=mirna_exp+runif(nrow(tmp),min=0,max=0.001)) ->tmp
  broom::tidy(cor.test(tmp$gene_exp,tmp$mirna_exp,method = c("spearman")),
              warning =function(e) 2 ,
              error=function(e) 1) -> tmp.spm
  if(length(tmp.spm)!=1){
    return(tmp.spm)
  }
}

# calculation 
net_miRNA_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=net_gene_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="name","gene"="symbol") -> net_miR_gene_spm_cor

net_miR_gene_spm_cor %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(mirna_ppar_network,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene.x==gene.y) %>%
  dplyr::select(regulator,gene.x,Cor,p.value) -> net_all_cor

net_all_cor %>%
  readr::write_tsv(file.path(FFL_data_path,"miRNA_PPAR_network","miRNA_PPAR_correlation.txt"))
# correlation filter  
net_all_cor %>%
  dplyr::mutate(label=ifelse(p.value<=0.05 & abs(Cor)>0.3,"*","")) %>%
  dplyr::mutate(`-log10(P)`=-log10(p.value)) %>%
  dplyr::mutate(`-log10(P)`=ifelse(`-log10(P)`=="Inf" |`-log10(P)`>5,5,`-log10(P)`)) %>%
  dplyr::rename("gene"="gene.x")-> ready_draw


### plot ----
ready_draw %>%
  dplyr::group_by(regulator) %>%
  dplyr::mutate(cor_sum=sum(Cor)) %>%
  dplyr::arrange(cor_sum) %>%
  dplyr::select(regulator,cor_sum) %>%
  dplyr::ungroup() %>%
  unique() -> y_rank
ready_draw %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(cor_sum=sum(Cor)) %>%
  dplyr::arrange(cor_sum) %>%
  dplyr::select(gene,cor_sum) %>%
  dplyr::ungroup() %>%
  unique() -> x_rank
library(ggplot2)
library(grid)
CPCOLS <- c("red", "white", "#1C86EE")

ready_draw %>%
  ggplot(aes(x=regulator,y=gene)) +
  geom_tile(aes(fill = Cor),color="white") +
  guides(fill = guide_colorbar(title.position = "left",
                               title.theme = element_text(angle = 90))) +
  scale_fill_gradient2(
    name = "Correlation", 
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-0.2,0,0.2)
  ) +
  scale_x_discrete(limit = y_rank$regulator) +
  scale_y_discrete(limit = x_rank$gene) +
  geom_text(aes(label=label)) +
  xlab("Upregulated miRNA") +
  ylab("TSGs targeted by CBX2 and EZH2") +
  coord_flip() +
  theme(#legend.position = "bottom",
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2),
    axis.text.x = element_text(size = 10,colour = "black",angle = 90,hjust = 1),
    axis.text.y = element_text(size = 10,colour = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black") ,
    plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(result_path,"Figure4/Figure5","PPAR_miRNA_Cor.pdf"),device = "pdf",width = 4,height = 4)
ggsave(file.path(result_path,"Figure4/Figure5","PPAR_miRNA_Cor.tiff"),device = "tiff",width = 4,height = 4)
