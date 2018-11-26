# TSG and oncogene statistic ----------------------------------------------
# E zhou -----
TSG_onco_data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"

# Hust ------
TSG_onco_data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"

TSG <- readr::read_tsv(file.path(TSG_onco_data_path,"TSG.source_clear(at least one evidence-no confuse).tsv")) %>%
  dplyr::mutate(hallmark="TSG")
oncogene <- readr::read_tsv(file.path(TSG_onco_data_path,"oncogene.source_clear(at least one evidence-no confuse).tsv")) %>%
  dplyr::mutate(hallmark="oncogene")
confuse_gene <- readr::read_tsv(file.path(TSG_onco_data_path,"confuse_gene.source_clear(at least one evidence).tsv")) %>%
  dplyr::mutate(hallmark="confuse") %>%
  dplyr::select(symbol,hallmark) %>%
  dplyr::rename("SYMBOL"="symbol")
TSG %>%
  rbind(oncogene) %>%
  dplyr::select(SYMBOL,hallmark) %>%
  rbind(confuse_gene) -> all_cancer_relate_genes_noconfused


# Regulation of TFs in EZH2/CBX2 targets to to TSGs and PPAR --------------------------------------------------------
### load miRNa regulation data -----
### >>>> construc FFL in server 1:/home/huff/LUAD_cancer/FFL_quantification_data/EZH2_CBX2_targets
# Ezhou ------
FFL_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets-180830"

# HUST -------
FFL_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets-180830"

net <- readr::read_tsv(file.path(FFL_path,"TF2gene"),col_names = F) %>%
  dplyr::rename("Sources"="X1","Targets"="X2")

### spearman correlation analysis filter ----

# data path config
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# load expression 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 

# get expression 
gene_exp %>%
  dplyr::filter(symbol %in% c(net$Targets)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_gene_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(net$Sources)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> net_TF_exp.gather

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
net_TF_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=net_gene_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="symbol","gene"="symbol1") -> net_TF_gene_spm_cor

net_net <- net
net_net %>%
  dplyr::rename("regulator"="Sources","gene"="Targets") -> net_net
net_TF_gene_spm_cor %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(net_net,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene.x==gene.y) %>%
  dplyr::select(regulator,gene.x,Cor,p.value) -> net_all_cor

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
  unique() -> tf_rank

readr::read_tsv(file.path(FFL_path,"gene.id")) %>%
  dplyr::mutate(n=c(1:26)) %>%
  dplyr::mutate(n=ifelse(n<=13,"TSG","PPAR")) %>%
  dplyr::rename("gene"="SYMBOL")-> gene.id
ready_draw %>%
  dplyr::inner_join(gene.id,by="gene") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(cor_sum=sum(Cor)) %>%
  dplyr::arrange(n,cor_sum) %>%
  dplyr::select(gene,cor_sum) %>%
  dplyr::ungroup() %>%
  unique() -> gene_rank
library(ggplot2)
library(grid)
CPCOLS <- c("red", "white", "#1C86EE")

ready_draw %>%
  ggplot(aes(y=regulator,x=gene)) +
  geom_tile(aes(fill = Cor),color="white") +
  guides(fill = guide_colorbar(title.position = "left",
                               title.theme = element_text(angle = 90))) +
  scale_fill_gradient2(
    name = "Correlation", 
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-0.1,0,0.3,0.6)
  ) +
  scale_x_discrete(limit = gene_rank$gene) +
  scale_y_discrete(limit = tf_rank$regulator) +
  geom_text(aes(label=label)) +
  xlab("PPAR signaling   TSGs") +
  ylab("TFs targeted by CBX2 and EZH2") +
  theme(#legend.position = "bottom",
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2),
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black") ,
    plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(result_path,"Figure4/Figure5","TSG-PPAR-targets_TF_Cor.pdf"),device = "pdf",width = 8,height = 3)
ggsave(file.path(result_path,"Figure4/Figure5","TSG-PPAR-targets_TF_Cor.tiff"),device = "tiff",width = 8,height = 3)

### add fold change info ----
ready_draw %>%
  dplyr::select(regulator,gene) -> network.cor.filter
data_path_3 <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/noiseq_no_cutoff_result"
TF_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_TF_cpm_1_noFDR")) 
progene_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_ProGene_cpm_1_noFDR"))
mirna_nofil <- readr::read_tsv(file.path(data_path_3,"NOISeq_DE_mirna_noFC_cpm30.mirnaid"))
rbind(TF_nofil,progene_nofil) %>%
  rbind(mirna_nofil) %>%
  dplyr::rename("SYMBOL"="gene_id") -> all_gene_nofil

net %>%
  tidyr::gather(key="type",value="Gene") %>%
  dplyr::mutate(type=ifelse(type=="Sources","1","2")) -> attribute

gene.id %>%
  dplyr::rename("SYMBOL"="gene","Class"="n") -> gene.id
attribute %>%
  dplyr::filter(Gene %in% c(net_all_cor$regulator,net_all_cor$gene.x)) %>%
  dplyr::rename("SYMBOL"="Gene") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,type,log2FC) %>%
  dplyr::left_join(all_cancer_relate_genes_noconfused,by="SYMBOL") %>%
  dplyr::left_join(gene.id,by="SYMBOL") -> attribute.fc

attribute.fc %>%
  readr::write_tsv(file.path(FFL_path,"attribute_fc_tsg_multuiple-source.txt"))
net_net %>%
  dplyr::mutate(regu_type=1) %>%
  dplyr::select(-X3) %>%
  readr::write_tsv(file.path(FFL_path,"network_multuiple-source.txt"))
