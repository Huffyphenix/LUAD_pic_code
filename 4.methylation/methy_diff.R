.libPaths("E:/library")
# load package ------------------------------------------------------------

library(magrittr)
library(ggplot2)
# data path ---------------------------------------------------------------
# E Zhou -----
data_path <- "H:/data/TCGA/LUAD_methylation"
gene_list_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# Hust -----
data_path <- "G:/data/TCGA/LUAD_methylation"
gene_list_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
miRNA_exp_path <- "G:/data/TCGA/TCGA_data"

# output data -------------------------------------------------------------
# E Zhou -----
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# Hust -----
out_path_sup <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# load data ---------------------------------------------------------------

luad_meth <- readr::read_tsv(file.path(data_path,"LUAD_paired_methylation.txt"))


gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 

genelist_pro <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
genelist_TF <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))


# data filter -------------------------------------------------------------
fn_rank <- function(data){
  data %>%
    dplyr::arrange(Genomic_Coordinate) %>%
    dplyr::mutate(Genomic_Direction=1:nrow(data)) %>%
    dplyr::mutate(Genomic_Direction=as.character(Genomic_Direction))  -> .out
  return(.out)
}

fn_median <- function(data){
  median(data$methy) -> .out
  return(.out)
}

luad_meth %>%
  dplyr::filter(Gene_Symbol %in% c(genelist_pro$gene_id)) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::mutate(group=ifelse(substr(sample,6,6)==1,"N","T")) -> genelist_pro.methy
luad_meth %>%
  dplyr::filter(Gene_Symbol %in%  c(genelist_pro$gene_id)) %>%
  dplyr::select(tag,Genomic_Coordinate,Gene_Symbol) %>%
  dplyr::mutate(Genomic_Coordinate=as.numeric(Genomic_Coordinate)) %>%
  tidyr::nest(-Gene_Symbol) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=purrr::map(data,fn_rank)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Gene_Symbol) -> genelist_pro.tag_posi
genelist_pro.methy %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group,-Gene_Symbol) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(genelist_pro.tag_posi,by="tag") %>%
  dplyr::ungroup() -> genelist_pro.median.methy
  
luad_meth %>%
  dplyr::filter(Gene_Symbol %in% c(genelist_TF$gene_id)) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::mutate(group=ifelse(substr(sample,6,6)==1,"N","T")) -> genelist_tf.methy
luad_meth %>%
  dplyr::filter(Gene_Symbol %in%  c(genelist_TF$gene_id)) %>%
  dplyr::select(tag,Genomic_Coordinate,Gene_Symbol) %>%
  dplyr::mutate(Genomic_Coordinate=as.numeric(Genomic_Coordinate)) %>%
  tidyr::nest(-Gene_Symbol) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=purrr::map(data,fn_rank)) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Genomic_Coordinate,-tag1,-Genomic_Coordinate1,-Gene_Symbol)-> genelist_tf.tag_posi
genelist_tf.methy %>%
  dplyr::select(-Gene_Symbol) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(genelist_tf.tag_posi,by="tag") -> genelist_tf.median.methy

# for CBX2 and EZH2 ------
data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/图"
cbx2 <- readr::read_tsv(file.path(data_path,"CBX2_promoter.txt"))
ezh2 <- readr::read_tsv(file.path(data_path,"EZH2_promoter.txt"))
gene_list <- c("EZH2","CBX2")
tags <- rbind(cbx2,ezh2)
# for other gene list ------
gene_list=c("PPARG")


pparg <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "PPARG") %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[10]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
FABP5 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "FABP5") %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::select(tag,Genomic_Coordinate) %>%
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[11] & Genomic_Coordinate>=Genomic_Coordinate[2]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
OLR1 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "OLR1") %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::select(tag,Genomic_Coordinate) %>%
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[4]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
rbind(pparg,FABP5,OLR1) -> tags


gene_list=c("CLDN11","RBMS3")
cldn11 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "CLDN11") %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[11]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
rbms3 <- luad_meth %>%
  dplyr::filter(Gene_Symbol %in% "RBMS3") %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::filter(Genomic_Coordinate<=Genomic_Coordinate[7]) %>% # according to gene strand(+,-),and the infomation of promoter from MEXPRESS database.
  dplyr::select(tag)
rbind(cldn11,rbms3) -> tags

# get promoter tags from mEXPRESS ------
luad_meth %>%
  dplyr::filter(tag %in% tags$tag) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::mutate(group=ifelse(substr(sample,6,6)==1,"N","T")) -> gene_list.methy
luad_meth %>%
  dplyr::filter(tag %in% tags$tag) %>%
  dplyr::arrange(Genomic_Coordinate) %>% 
  dplyr::select(tag,Genomic_Coordinate,Gene_Symbol) %>%
  dplyr::mutate(Genomic_Coordinate=as.numeric(Genomic_Coordinate)) %>%
  tidyr::nest(-Gene_Symbol) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=purrr::map(data,fn_rank)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=paste("P",Genomic_Direction,sep="_"))-> gene_list.tag_posi
gene_list.methy %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group,-Gene_Symbol) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::ungroup()-> gene_list.median.methy

# calculation -------------------------------------------------------------

# EZH2 and CBX2 -----------------------------------------------------------
gene_list.methy %>%
  tidyr::drop_na(methy) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  dplyr::group_by(tag) %>%
  dplyr::do(
    broom::tidy(
      wilcox.test(methy ~ group, data = .)
    )
  ) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.01,"**",sig)) %>%
  dplyr::inner_join(gene_list.methy,by="tag") %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::mutate(laby=max(as.numeric(methy))+0.05) %>%
  dplyr::select(tag,Gene_Symbol,sig,Genomic_Direction,laby) %>%
  dplyr::ungroup() %>%
  unique()-> gene_list.ttest

gene_list.methy %>%
  dplyr::inner_join(gene_list.tag_posi,by="tag") %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::drop_na(methy) %>%
  ggplot(aes(x=Genomic_Direction,y=methy,color=group)) +
  geom_boxplot() +
  geom_point(position = "jitter",size=1) +
  geom_text(data=gene_list.ttest,mapping=aes(x=Genomic_Direction,y=laby,label=sig),color="black") +
  scale_color_manual(values = c("#5CACEE", "#FF82AB")) +
  # scale_x_discrete(limit = EZH2_CBX2.tag_posi$Genomic_Direction[1:27]) +
  facet_wrap(~Gene_Symbol,scales = "free") +
  ylab(latex2exp::TeX("Methylation ($\\beta$ value)")) +
  xlab("Promoter") +
  theme(
    axis.line = element_line(color = "black"),
    panel.background  = element_rect(fill = "white", color = "black"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    # axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20),
    axis.ticks.x = element_blank(),
    legend.position = c(0.9,0.8)
  ) -> EZH2_CBX2.box;EZH2_CBX2.box
ggsave(file.path(out_path_fig,"Figure3","Figure3B.Methy_histone.pdf"),EZH2_CBX2.box,device = "pdf",width = 6,height = 3)
ggsave(file.path(out_path_fig,"Figure3","Figure3B.Methy_histone.tiff"),EZH2_CBX2.box,device = "tiff",width = 6,height = 3)
ggsave(file.path(out_path_fig,"Figure4/Figure5","PPAR_Methy_boxplot.pdf"),EZH2_CBX2.box,device = "pdf",width = 10,height = 6)
ggsave(file.path(out_path_fig,"Figure4/Figure5","PPAR_Methy_boxplot.tiff"),EZH2_CBX2.box,device = "tiff",width = 10,height = 6)
ggsave(file.path(out_path_fig,"Figure4/Figure5","CLDN11_RBMS3_Methy_boxplot.pdf"),EZH2_CBX2.box,device = "pdf",width = 5,height = 3)
ggsave(file.path(out_path_fig,"Figure4/Figure5","CLDN11_RBMS3_Methy_boxplot.tiff"),EZH2_CBX2.box,device = "tiff",width = 5,height = 3)

EZH2_CBX2.median.methy %>%
  dplyr::mutate(Genomic_Direction=as.numeric(Genomic_Direction)) %>%
  unique() %>%
  ggplot(aes(x=Genomic_Direction,y=median,color=group)) +
  geom_line() +
  scale_x_discrete(limit = EZH2_CBX2.tag_posi$Genomic_Direction[1:27]) +
  facet_grid(~Gene_Symbol) +
  theme(
    axis.line = element_line(color = "black"),
    panel.background  = element_rect(fill = "white", color = "grey"),
    panel.grid = element_line(colour = "grey"),
    # axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    text = element_text(size = 20)
    ) -> EZH2_CBX2.line;EZH2_CBX2.line
ggsave(file.path(out_path_fig,"Figure3","Figure3B.Methy_histone_line.pdf"),EZH2_CBX2.line,device = "pdf",width = 10,height = 6)


# targets -----------------------------------------------------------------

genelist_pro.methy %>%
  tidyr::drop_na(methy) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  # dplyr::filter(tag=="cg00008737") %>%
  dplyr::mutate(sample=substr(sample,1,4)) %>%
  dplyr::group_by(tag,Gene_Symbol,sample,group) %>%
  dplyr::mutate(methy_mean=mean(methy)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-methy) %>%
  unique() %>%
  tidyr::spread(group,methy_mean) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::do(
    broom::tidy(
      t.test(.$`T`, .$N, data = ., paired=TRUE)
    )
  ) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.01,"**",sig)) %>%
  dplyr::mutate(`diff(T-N)` = estimate) %>%
  dplyr::rename("ttestPvalue"="p.value","ttestMethod"="method","diffSig"="sig") -> genelist_pro.paired_ttest

genelist_pro.paired_ttest %>%
  dplyr::filter(`diff(T-N)`<0 & p.value<=0.05) -> genelist_diff_in_methyTN


gene_exp %>%
  dplyr::filter(symbol %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol,.key="exp") %>%
  dplyr::rename("Gene_Symbol"="symbol")-> gene_exp.gather
luad_meth %>%
  dplyr::filter(Gene_Symbol %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  dplyr::group_by(Gene_Symbol,sample) %>%
  dplyr::mutate(methy_mean=mean(methy)) %>%
  dplyr::select(-tag,-methy) %>%
  unique() %>%
  dplyr::ungroup() %>% 
  tidyr::nest(-Gene_Symbol,.key="methy") -> genelist_methy.gather

fn_spm <- function(.exp,.methy){
  .exp %>%
    dplyr::inner_join(.methy,by="sample") %>%
    dplyr::mutate(methy_mean=methy_mean+runif(n(),min=0,max=0.001)) %>%
    dplyr::mutate(gene_exp=gene_exp+runif(n(),min=0,max=0.001)) -> .tmp
  broom::tidy(
      cor.test(.tmp$methy_mean,.tmp$gene_exp,method = "spearman")
    ) -> .out
  return(.out)
}
gene_exp.gather %>%
  dplyr::inner_join(genelist_methy.gather,by="Gene_Symbol") %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(spm=purrr::map2(exp,methy,fn_spm)) %>%
  dplyr::select(-(exp:methy)) %>%
  tidyr::unnest() %>%
  dplyr::rename("spmCor"="estimate") %>%
  dplyr::rename("spm_pvalue"="p.value","CorMethod"="method") -> genelist_exp_methy_cor
genelist_exp_methy_cor %>%
  dplyr::filter(spm_pvalue<=0.05) %>%
  dplyr::arrange(spmCor) %>% .$Gene_Symbol -> genesymbol_rank

genelist_exp_methy_cor %>%
  dplyr::inner_join(genelist_pro.paired_ttest,by="Gene_Symbol") %>%
  dplyr::select(Gene_Symbol,spmCor,spm_pvalue,ttestPvalue,`diff(T-N)`) -> genelist_spm_diff

CPCOLS <- c("red", "white", "blue")
genelist_exp_methy_cor %>%
  dplyr::filter(spm_pvalue<=0.05) %>%
  dplyr::mutate(spm_pvalue=ifelse(spm_pvalue==0,0.0000001,spm_pvalue)) %>%
  dplyr::mutate(log10pvalue=-log10(spm_pvalue)) %>%
  dplyr::mutate(x="Spearman Cor") %>%
  ggplot(aes(x=x,y=Gene_Symbol)) +
  geom_point(aes(size=log10pvalue,color=spmCor)) +
  scale_y_discrete(limit = genesymbol_rank) +
  scale_color_gradient2(
    name = "Spearman Cor", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1]
  ) +
  scale_size_continuous(
    name = "-log10(Pvalue)"
  ) +
  theme(legend.position = "bottom",
       panel.background = element_rect(colour = "black", fill = "white"),
       panel.grid = element_line(colour = "grey", linetype = "dashed"),
       panel.grid.major = element_line(
         colour = "grey",
         linetype = "dashed",
         size = 0.2),
       axis.text.x = element_text(size = 10),
       axis.text.y = element_text(size = 10),
       legend.text = element_text(size = 10),
       legend.title = element_text(size = 12),
       legend.key = element_rect(fill = "white", colour = "black") ,
       plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(out_path_fig,"Figure4","Figure4B.methy_diff.pdf"),device = "pdf",width = 8,height = 8)

genelist_pro.paired_ttest %>%
  dplyr::filter(Gene_Symbol %in% genesymbol_rank) %>%
  dplyr::mutate(ttestPvalue=ifelse(ttestPvalue==0,0.0000001,ttestPvalue)) %>%
  dplyr::mutate(log10pvalue=-log10(ttestPvalue)) %>%
  dplyr::mutate(x="Methy_Diff(T-N)") %>%
  ggplot(aes(x=x,y=Gene_Symbol)) +
  geom_point(aes(size=log10pvalue,color=`diff(T-N)`)) +
  scale_y_discrete(limit = genesymbol_rank) +
  scale_color_gradient2(
    name = "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1]
  ) +
  scale_size_continuous(
    name = "-log10(Pvalue)"
  ) +
  theme(legend.position = "bottom",
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = "white", colour = "black") ,
        plot.title = element_text(size=20)
  ) -> p;p
ggsave(file.path(out_path_fig,"Figure4","Figure4B.methy_Cor.pdf"),device = "pdf",width = 8,height = 8)
