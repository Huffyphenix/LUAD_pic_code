# 上调的miRNA + 下调的TF +CBX6 、7构建FFL
# .3 /project/huff/huff/TCGA_lungcancer_data/pathway_enrichment_FFL/CBX_FFL

# 添加差异表达和TSG等信息
FFL_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX4678/CBX_FFL"
# TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"

# .libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")

# load data ---------------------------------------------------------------
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

attribute <- readr::read_tsv(file.path(FFL_data_path,"attribute.txt"),col_names = F)
network <- readr::read_tsv(file.path(FFL_data_path,"network.txt"),col_names = F) 

network %>% 
  dplyr::rename("source"="X1","target"="X2","regulate_type"="X3") %>%
  readr::write_tsv(file.path(FFL_data_path,"network.txt"))

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

# add hallmark info -------------------------------------------------------

attribute %>%
  dplyr::rename("SYMBOL"="X1","gene_type"="X2") %>%
  dplyr::left_join(all_cancer_relate_genes_noconfused,by="SYMBOL") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::mutate(mark=ifelse(hallmark=="TSGene",2,1)) %>%
  dplyr::mutate(mark=ifelse(is.na(hallmark),0,mark)) %>%
  dplyr::select(SYMBOL,gene_type,log2FC,mark) -> attribute.hallmark
attribute.hallmark %>%
  readr::write_tsv(file.path(FFL_data_path,"attribute.multi-source.hallmark-added.txt"))

# cytoscape构建网络
# filter ---
# 删除与CBX6，7无直接调控作用的调控因子
# 只保留有直接作用的
# 删除与CBX表达量变化趋势不一致的调控因子，即CBX表达量上调，则只保留调控它的上调的TF和下调的miRNA
regulatory <- readr::read_tsv(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/CBX4678/CBX_FFL","cytoscape_filter_net.txt"),col_names = F) %>%
  dplyr::filter(X2 %in% c("CBX2","CBX4","CBX6","CBX7","CBX8"))

# correlation analysis
# data path config
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# load expression 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
mirna_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 

gene_exp %>%
  dplyr::filter(symbol %in% c(regulatory$X2)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_gene_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(regulatory$X1)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_TF_exp.gather
mirna_exp %>%
  dplyr::filter(name %in% c(regulatory$X1)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> net_mirna_exp.gather

# correlation for gene2gene
fn_get_spm_a_g <- function(m_data,g_exp){
  # print(m_data)
  g_exp %>%
    dplyr::group_by(symbol) %>%
    # dplyr::filter(symbol %in% c("CBX2")) %>%
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
  set.seed(10000)
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
# function for mirna2gene
fn_get_spm_a <- function(m_data,g_exp){
  # print(m_data)
  g_exp %>%
    dplyr::group_by(symbol) %>%
    # dplyr::filter(symbol %in% c("CBX2")) %>%
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
  set.seed(10000)
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

# correlation calculation  -----
net_mirna_exp.gather %>%
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

regulatory %>%
  dplyr::rename("regulator"="X1","gene"="X2") -> net_net
net_miR_gene_spm_cor %>%
  rbind(net_TF_gene_spm_cor) %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(net_net,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene.x==gene.y) %>%
  dplyr::select(regulator,gene.x,Cor,p.value) -> net_all_cor

net_all_cor %>%
  readr::write_tsv(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/CBX4678/CBX_FFL","net_corrlation_analysis.tsv"))

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
  geom_tile(aes(fill=Cor),colour = "grey") +
  # geom_point(size=3) +
  guides(color=guide_colorbar(title.position="left")) +
  scale_x_discrete(limits = x_rank$gene) +
  scale_y_discrete(limits = y_rank$regulator) +
  ylab("Regulators") +
  # xlab("CBX2 and EZH2 targets (downreulate)") +
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
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.position = c(0.8,0.5),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20)
  ) -> p;p
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX4678/CBX_FFL/"
ggsave(file.path(out_path_fig,"network.correlation.pdf"),device = "pdf",width = 4,height = 3)
ggsave(file.path(out_path_fig,"network.correlation.tiff"),device = "tiff",width = 4,height = 3)
