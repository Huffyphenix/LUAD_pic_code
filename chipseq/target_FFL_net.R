# data path ---------------------------------------------------------------

genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
data_path<- "H:/data"
chip_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"

library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
# data manage -------------------------------------------------------------
genelist <- readr::read_tsv(file.path(genelist_path,"all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(Class=="Down")

Animal_TF <-  readr::read_tsv(file.path(data_path,"AnimalTFDB","Homo_sapiens_transcription_factors_gene_list.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls")) %>%
  .$Symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Animal_TF %>%
  dplyr::filter(! Entrez_ID %in% enzyme_lsit$ENTREZID) -> Animal_TF
Animal_TF$Entrez_ID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db) -> Animal_TF_entrez

genelist %>%
  dplyr::filter(entrez_id %in% Animal_TF_entrez$ENTREZID) -> DOWN_commone_targets.TF 
DOWN_commone_targets.TF %>%
  readr::write_tsv(file.path(genelist_path,"FFL/input","DOWN_commone_targets.TF"))
genelist %>%
  dplyr::filter(! entrez_id %in% Animal_TF_entrez$ENTREZID) -> DOWN_commone_targets.pro
DOWN_commone_targets.pro %>%
  readr::write_tsv(file.path(genelist_path,"FFL/input","DOWN_commone_targets.pro"))

# TSG and oncogene statistic ----------------------------------------------
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

genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(gene_id.x %in% TSGene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(gene_id.x %in% oncogene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC>0) %>%
  dplyr::filter(gene_id.x %in% oncogene_noconfuse$symbol) %>% nrow()
genelist %>%
  dplyr::filter(log2FC>0) %>%
  dplyr::filter(gene_id.x %in% TSGene_noconfuse$symbol) %>% nrow()


### >>>> construc FFL in server 1:/home/huff/LUAD_cancer/FFL_quantification_data/EZH2_CBX2_targets/
ffl_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets"
TF2miRNA <- readr::read_tsv(file.path(ffl_path,"TF2miRNA"),col_names = F)
TF2gene <- readr::read_tsv(file.path(ffl_path,"TF2gene"),col_names = F)
miRNA2TF <- readr::read_tsv(file.path(ffl_path,"miRNA2TF"),col_names = F)
miRNA2gene <- readr::read_tsv(file.path(ffl_path,"miRNA2gene"),col_names = F)
attribute <- readr::read_tsv(file.path(ffl_path,"attribute.txt"),col_names = F)

# filter ------------------------------------------------------------------

# only remain E2F1 and SOX4 as upstream of miRNAs, cause other TFs are all downregulate because of EZH2 and CBX2.
# TF2miRNA %>%
#   dplyr::filter(X1 %in% c("E2F1","SOX4")) %>%
#   dplyr::mutate(X4=1) -> E2F1_SOX4_2_miRNA


miRNA2TF %>%
  dplyr::filter(X2 %in% unique(DOWN_commone_targets.TF$gene_id.x)) %>%
  dplyr::mutate(X3="notshow") %>%
  dplyr::mutate(X4=3) -> miRNA2TF.filter

miRNA2gene %>%
  # dplyr::filter(X1 %in% unique(E2F1_SOX4_2_miRNA$X2)) %>%
  dplyr::mutate(X3="notshow") %>%
  dplyr::mutate(X4=4) -> miRNA2gene.filter

TF2gene %>%
  dplyr::filter(! X1 %in% c("E2F1","SOX4")) %>%
  dplyr::filter(X2 %in% unique(miRNA2gene.filter$X2)) %>%
  dplyr::mutate(X4=2)-> TF2gene.filter

# combine -----------------------------------------------------------------

miRNA2TF.filter %>%
  rbind(miRNA2gene.filter) %>%
  rbind(TF2gene.filter) -> network

data.frame(gene=miRNA2TF.filter$X2,type=1) %>%
  rbind(data.frame(gene=miRNA2TF.filter$X1,type=2)) %>%
  rbind(data.frame(gene=miRNA2gene.filter$X1,type=2)) %>%
  rbind(data.frame(gene=miRNA2gene.filter$X2,type=3)) %>%
  rbind(data.frame(gene=TF2gene.filter$X1,type=1)) %>%
  rbind(data.frame(gene=TF2gene.filter$X2,type=3)) %>%
  dplyr::mutate(gene=as.character(gene)) %>%
  unique()  -> attribute

# spearman correlation analysis filter ------------------------------------

# data path config
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# load expression 
miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 

# data manage
net_TF <- attribute %>% dplyr::filter(type==1)
net_miR <- attribute %>% dplyr::filter(type==2)
net_gene <- attribute %>% dplyr::filter(type==3)

miRNA_exp %>%
  dplyr::filter(name %in% c(net_miR$gene)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> net_miRNA_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(net_TF$gene)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_TF_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(net_gene$gene,net_TF$gene)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_gene_exp.gather

# function for mirna2gene
fn_get_spm_a <- function(m_data,g_exp){
  # print(m_data)
  g_exp %>%
    dplyr::group_by(symbol) %>%
    # dplyr::filter(gene_id %in% c("RBMS3")) %>%
    dplyr::mutate(spm=purrr::map(data,m_data=m_data,fn_get_spm_b)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::ungroup() %>%
    dplyr::rename(Cor=estimate)-> .out
  return(.out)
}

fn_get_spm_b <- function(g_data,m_data){
  # print(g_data)
  # print(m_data)
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

# correlation for gene2gene
fn_get_spm_a_g <- function(m_data,g_exp){
  # print(m_data)
  g_exp %>%
    dplyr::group_by(symbol) %>%
    # dplyr::filter(gene_id %in% c("RBMS3")) %>%
    dplyr::mutate(spm=purrr::map(data,m_data=m_data,fn_get_spm_b_g)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::ungroup() %>%
    dplyr::rename(Cor=estimate)-> .out
  return(.out)
}

fn_get_spm_b_g <- function(g_data,m_data){
  # print(g_data)
  # print(m_data)
  g_data %>%
    dplyr::inner_join(m_data,by="sample") ->tmp
  tmp %>%
    dplyr::mutate(gene_exp.x=gene_exp.x+runif(nrow(tmp),min=0,max=0.001)) %>%
    dplyr::mutate(gene_exp.y=gene_exp.y+runif(nrow(tmp),min=0,max=0.001)) ->tmp
  broom::tidy(cor.test(tmp$gene_exp.x,tmp$gene_exp.y,method = c("spearman")),
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

net_TF_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=net_gene_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="symbol","gene"="symbol1") -> net_TF_gene_spm_cor

net_net <- network
net_net %>%
  dplyr::rename("regulator"="X1","gene"="X2") -> net_net
net_miR_gene_spm_cor %>%
  rbind(net_TF_gene_spm_cor) %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(net_net,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene.x==gene.y) %>%
  dplyr::select(regulator,gene.x,Cor,p.value,X3,X4) -> net_all_cor

# correlation filter  
net_all_cor %>%
  dplyr::filter(p.value<=0.05 & abs(Cor)>0.3) %>%
  dplyr::mutate(`-log10(P)`=-log10(p.value)) %>%
  dplyr::mutate(`-log10(P)`=ifelse(`-log10(P)`=="Inf" |`-log10(P)`>5,5,`-log10(P)`)) %>%
  dplyr::rename("gene"="gene.x")-> ready_draw

# draw picture
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
CPCOLS <- c("red", "white", "blue")
ready_draw %>%
  ggplot(aes(x=gene,y=regulator,color=Cor)) +
  geom_tile(fill=c("#EDEDED"),colour = "grey") +
  geom_point(size=3) +
  guides(color=guide_colorbar(title.position="left")) +
  scale_x_discrete(limits = x_rank$gene) +
  scale_y_discrete(limits = y_rank$regulator) +
  ylab("Regulators") +
  xlab("CBX2 and EZH2 targets (downreulate)") +
  scale_color_gradient2(
    name = "Correlation", # "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1],
    breaks=c(-0.3,0,0.3,0.6)
  ) +
  theme(
    # legend.position = "bottom",
    panel.background = element_rect(colour = "black", fill = "white"),
    # panel.grid = element_line(colour = "grey", linetype = "dashed"),
    # panel.grid.major = element_line(
    #   colour = "grey",
    #   linetype = "dashed",
    #   size = 0.2
    # ),
    
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, angle = 40, size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.position = c(0.8,0.5),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20)
  ) -> p;p
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"
ggsave(file.path(out_path_fig,"Figure4/Figure5","Figue S7C.network.correlation.pdf"),device = "pdf",width = 10,height = 4)
ggsave(file.path(out_path_fig,"Figure4/Figure5","Figue S7C.network.correlation.tiff"),device = "tiff",width = 10,height = 4)

# add fold change info ----------------------------------------------------
ready_draw %>%
  dplyr::select(regulator,gene,X3,X4) -> network.cor.filter
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute %>%
  dplyr::filter(gene %in% c(network.cor.filter$regulator,network.cor.filter$gene)) %>%
  dplyr::rename("SYMBOL"="gene") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,type,log2FC) -> attribute.fc

# add TSG info ------------------------------------------------------------

attribute.fc %>%
  dplyr::left_join(all_cancer_relate_genes_noconfused,by="SYMBOL") %>%
  dplyr::mutate(mark=ifelse(hallmark=="TSGene",2,1)) %>%
  dplyr::mutate(mark=ifelse(is.na(hallmark),0,mark)) %>%
  dplyr::rename("symbol"="SYMBOL") %>%
  dplyr::select(symbol,type,log2FC,mark) -> attribute.fc.tsg
  

attribute.fc.tsg %>%
  readr::write_tsv(file.path(genelist_path,"FFL/EZH2_CBX2_targets/√network.filter","attribute_fc_tsg_multuiple-source.txt"))
network.cor.filter %>%
  readr::write_tsv(file.path(genelist_path,"FFL/EZH2_CBX2_targets/√network.filter","network_multuiple-source.txt"))



### mirna network select hub nodes from cytoNCA------
cytonca <- readr::read_tsv(file.path(ffl_path,"√network.filter","cytoCNA.txt"),col_names=F)

cytonca %>%
  tidyr::separate(X3,c("group","Subgragh"),sep=": ") %>%
  tidyr::separate(X4,c("group1","Degree"),sep=": ") %>%
  tidyr::separate(X5,c("group2","Eigenvector"),sep=": ") %>%
  tidyr::separate(X6,c("group3","Information"),sep=": ") %>%
  tidyr::separate(X7,c("group4","LAC"),sep=": ") %>%
  tidyr::separate(X8,c("group5","Betweenness"),sep=": ") %>%
  tidyr::separate(X9,c("group6","Closeness"),sep=": ") %>%
  tidyr::separate(X10,c("group7","Network"),sep=": ") %>%
  dplyr::select(X2, Subgragh,Degree,Eigenvector,Information,LAC,Betweenness,Closeness,Network) %>%
  tidyr::gather(-X2,key="group",value="value") %>%
  dplyr::mutate(value=as.numeric(value)) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(rank=rank(value)) %>% # use average rank when a ties happened
  dplyr::group_by(X2) %>%
  dplyr::mutate(rank_g=sum(rank)) %>%
  dplyr::ungroup() -> cytonca_rank


cytonca_rank %>%
  dplyr::select(-rank_g,-rank) %>%
  tidyr::spread(key=group,value=value) -> cytonca_value

cytonca_rank %>%
  dplyr::select(-value) %>%
  dplyr::mutate(group=paste(group,"rank",sep="_")) %>%
  tidyr::spread(key=group,value=rank) %>%
  dplyr::inner_join(cytonca_value,by="X2") %>%
  dplyr::rename("Gene"="X2","Rank_sum"="rank_g")-> cytonca_rank_all

cytonca_rank_all %>%
  readr::write_tsv(file.path(ffl_path,"√network.filter","cytoCNA_rank.txt"))

# genes out of control of CNV, methylation, and miRNA regulation -----
# 1 methylation
mthy_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"
methy_regu_gene <- readr::read_tsv(file.path(mthy_path,"Figure4/Figure5","genes_regulate_by_methy.tsv"))

# 2 CNV 
Down_common_targets.CNV <- readr::read_tsv(file.path(mthy_path,"Figure4/Figure5","Down_common_targets.CNV-percent.tsv"))
Down_common_targets.CNV %>%
  dplyr::filter(type== "Deletion") %>%
  dplyr::filter(Percent > 5) -> Down_common_targets.CNV_dele_5

# 3 miRNA regulation
attribute.fc.tsg %>%
  dplyr::filter(type %in% c(1,3)) %>%
  dplyr::filter(! symbol %in% c("E2F1","SOX4")) -> miRNA_regu_genes
# 
# cytonca_rank_all %>%
#   dplyr::arrange(desc(Rank_sum)) %>%
#   dplyr::filter(Rank_sum>Rank_sum[68]) -> miRNA_regu_genes

# 4 filter
genelist %>%
  dplyr::filter(log2FC<0) %>%
  dplyr::filter(! gene_id.x %in% methy_regu_gene$symbol) %>%
  dplyr::filter(! gene_id.x %in% Down_common_targets.CNV_dele_5$SAMPLE_ID) %>%
  dplyr::filter(! gene_id.x %in% miRNA_regu_genes$symbol) %>%
  dplyr::filter(gene_id.x %in% TSG$SYMBOL)
