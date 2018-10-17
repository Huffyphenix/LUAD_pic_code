# TSG and oncogene statistic ----------------------------------------------
# E zhou ----
TSG_onco_data_path <- "F:/?业募?????/ENCODE-TCGA-LUAD/TS and oncogene source"

# HUST ----
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


# EZH2 CBX2 downregulated targets -----------------------------------------
# E zhou ----
chip_path <- "F:/?业募?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"

# HUST ----
chip_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"

all_EHZ2_CBX2_common_targets.DE_info <- readr::read_tsv(file.path(chip_path,"common-targets-180426-new","all_EHZ2_CBX2_common_targets.DE_info"))

all_EHZ2_CBX2_common_targets.DE_info %>%
  dplyr::filter(Class=="Down") -> all_EHZ2_CBX2_common_targets.Down


# TSGenes in targets ------------------------------------------------------
TSG %>%
  dplyr::filter(SYMBOL %in% all_EHZ2_CBX2_common_targets.Down$gene_id.x) -> EZH2_CBX2_targets_TSG

# E zhou ----
result_path <- "F:/?业募?????/ENCODE-TCGA-LUAD/Figure/"

# HUST ----
result_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

EZH2_CBX2_targets_TSG %>%
  readr::write_tsv(file.path(result_path,"Figure4/Figure5","EZH2-CBX2_down_targets_TSG.tsv"))

# CNV analysis ------------------------------------------------------------
### load cnv data ----
data_path<- "H:/data/GSCALite/TCGA/cnv"
data_path<- "G:/data/GSCALite/TCGA/cnv"

luad_cnv <- readr::read_rds(file.path(data_path,"pancan34_cnv_percent.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()

### data filter ----
luad_cnv %>%
  dplyr::filter(symbol %in% EZH2_CBX2_targets_TSG$SYMBOL) %>%
  dplyr::select(symbol,a_homo,d_homo) %>%
  tidyr::gather(-symbol,key="type",value="CNV") %>%
  dplyr::mutate(Percent=CNV*100) %>%
  dplyr::mutate(type=ifelse(type=="a_homo","Amplification","Deletion"))-> cnv_plot_ready

### plot ----
cnv_plot_ready %>%
  dplyr::mutate(CNV=ifelse(type=="Deletion",-CNV,CNV)) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(CNV_sum=sum(CNV)) %>%
  dplyr::select(symbol,CNV_sum) %>%
  unique() %>%
  dplyr::arrange(CNV_sum) %>%
  dplyr::ungroup() ->symbol_rank

cnv_plot_ready %>%
  ggplot(aes(x=symbol,y=Percent,fill=type)) +
  geom_col(position = "stack",width = 0.3) +
  guides(fill=guide_legend(title = "")) +
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  scale_x_discrete(limits = symbol_rank$symbol) +
  ylab("CNV Percent (%)") +
  xlab("TSGs targeted by CBX2 and EZH2") +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    legend.position = c(0.7,0.3),
    legend.background = element_blank(),
    # axis.title.x = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(colour = "black",size = 12)
    # axis.text.x = element_text(angle = 30,hjust = 1)
  ) ->p;p
ggsave(file.path(result_path,"Figure4/Figure5","Down_TSG_targets.CNV-percent.pdf"),width = 3,height = 4)
ggsave(file.path(result_path,"Figure4/Figure5","Down_TSG_targets.CNV-percent.tiff"),width = 3,height = 4)

# DNA methylation ---------------------------------------------------------
### load methy data ----
# E zhou path
methy_data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/"

# HUST path 
methy_data_path <- "G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/"

methy <- readr::read_rds(file.path(methy_data_path,"pan33_allgene_methy_diff.simplification.rds.gz"))
methy_cor <- readr::read_rds(file.path(methy_data_path,"pancan34_all_gene_exp-cor-meth.rds.gz"))

### genelist filter ----
methy %>%
  dplyr::filter(cancer_types == "LUAD") %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% EZH2_CBX2_targets_TSG$SYMBOL) -> targets_methy_diff

methy_cor %>%
  dplyr::filter(cancer_types == "LUAD") %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% EZH2_CBX2_targets_TSG$SYMBOL) -> targets_methy_cor

targets_methy_diff %>%
  dplyr::full_join(targets_methy_cor,by="symbol") %>%
  dplyr::select(symbol,diff,spm) %>%
  tidyr::gather(-symbol,key="group",value="value") %>%
  dplyr::mutate(value=ifelse(is.na(value),0,signif(value,2))) %>%
  dplyr::mutate(label=ifelse(value==0,"NS",value)) %>%
  dplyr::mutate(group = ifelse(group=="diff","Diff. (T - N)","Cor."))-> targets_methy  # set NA value as 0

### plot ----
targets_methy %>%
  dplyr::filter(group == "spm") %>%
  dplyr::arrange(value) %>% .$symbol -> cor_rank.genesymbol

targets_methy_cor %>%
  dplyr::rename("value"="spm") %>%
  dplyr::mutate(group="Cor.") -> targets_methy_cor.pic

targets_methy_diff %>%
  dplyr::rename("value"="diff","logfdr" = "fdr") %>%
  dplyr::mutate(group="Diff. (T - N)") %>%
  dplyr::select(-direction) -> targets_methy_diff.pic

library(ggplot2)
library(grid)
CPCOLS <- c("red", "white", "#1C86EE")

targets_methy %>%
  ggplot(aes(x=group,y=symbol)) +
  geom_tile(aes(fill = value),color="white") +
  scale_y_discrete(limit = cor_rank.genesymbol) +
  guides(size = guide_legend(title.position = "left",
                             title.theme = element_text(angle = 90))) +
  scale_fill_gradient2(
    name = "Diff./Cor.", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
  ) +
  geom_text(aes(label=label)) +
  ylab("Symbol") +
  theme(#legend.position = "bottom",
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black") ,
    plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(result_path,"Figure4/Figure5","TSG_targets_methy_Cor-diff-gsca.pdf"),device = "pdf",width = 4,height = 6)
ggsave(file.path(result_path,"Figure4/Figure5","TSG_targets_methy_Cor-diff-gsca.tiff"),device = "tiff",width = 4,height = 6)

# miRNA regulation --------------------------------------------------------
### load miRNa regulation data -----
### >>>> construc FFL in server 1:/home/huff/LUAD_cancer/FFL_quantification_data/EZH2_CBX2_TSG_targets/

mirna_regulate_path <- "F:/?业募?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_TSG_targets"
net <- readr::read_tsv(file.path(mirna_regulate_path,"miRNA2gene"),col_names = F) %>%
  dplyr::rename("Sources"="X1","Targets"="X2")

### spearman correlation analysis filter ----

# data path config
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# load expression 
miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
# get expression 
miRNA_exp %>%
  dplyr::filter(name %in% c(net$Sources)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> net_miRNA_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(net$Targets)) %>%
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

net_net <- net
net_net %>%
  dplyr::rename("regulator"="Sources","gene"="Targets") -> net_net
net_miR_gene_spm_cor %>%
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
    breaks = c(-0.3,0,0.15)
  ) +
  scale_x_discrete(limit = y_rank$regulator) +
  scale_y_discrete(limit = x_rank$gene) +
  geom_text(aes(label=label)) +
  xlab("Upregulated miRNA") +
  ylab("TSGs targeted by CBX2 and EZH2") +
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
ggsave(file.path(result_path,"Figure4/Figure5","TSGtargets_miRNA_Cor.pdf"),device = "pdf",width = 8,height = 4)
ggsave(file.path(result_path,"Figure4/Figure5","TSGtargets_miRNA_Cor.tiff"),device = "tiff",width = 8,height = 4)

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

attribute %>%
  dplyr::filter(Gene %in% c(net_all_cor$regulator,net_all_cor$gene.x)) %>%
  dplyr::rename("SYMBOL"="Gene") %>%
  dplyr::left_join(all_gene_nofil,by="SYMBOL") %>%
  dplyr::select(SYMBOL,type,log2FC) %>%
  dplyr::left_join(EZH2_CBX2_targets_TSG,by="SYMBOL") %>%
  dplyr::select(-ENTREZID,-source) -> attribute.fc

attribute.fc %>%
  readr::write_tsv(file.path(mirna_regulate_path,"attribute_fc_tsg_multuiple-source.txt"))
net_net %>%
  dplyr::mutate(regu_type=1) %>%
  readr::write_tsv(file.path(mirna_regulate_path,"network_multuiple-source.txt"))
