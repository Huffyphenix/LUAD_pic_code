
# load package ------------------------------------------------------------

library(magrittr)
library(ggplot2)
# data path ---------------------------------------------------------------
data_path <- "H:/data/TCGA/LUAD_methylation"
gene_list_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"

# output data -------------------------------------------------------------
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# load data ---------------------------------------------------------------

luad_meth <- readr::read_tsv(file.path(data_path,"LUAD_paired_methylation.txt"))
luad_meth %>%
  tidyr::separate(gene,c("tag","symbol"),"_") -> luad_meth

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
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Genomic_Coordinate,-tag1,-Genomic_Coordinate1)-> genelist_pro.tag_posi
genelist_pro.methy %>%
  dplyr::select(-Gene_Symbol) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(genelist_pro.tag_posi,by="tag") -> genelist_pro.median.methy
  
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

luad_meth %>%
  dplyr::filter(Gene_Symbol %in% c("EZH2","CBX2")) %>%
  dplyr::select(-Genomic_Coordinate) %>%
  tidyr::gather(-tag,-Gene_Symbol,key="sample",value="methy") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::mutate(group=ifelse(substr(sample,6,6)==1,"N","T")) -> EZH2_CBX2.methy
luad_meth %>%
  dplyr::filter(Gene_Symbol %in% c("EZH2","CBX2")) %>%
  dplyr::select(tag,Genomic_Coordinate,Gene_Symbol) %>%
  dplyr::mutate(Genomic_Coordinate=as.numeric(Genomic_Coordinate)) %>%
  tidyr::nest(-Gene_Symbol) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(Genomic_Direction=purrr::map(data,fn_rank)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  dplyr::select(-Gene_Symbol) -> EZH2_CBX2.tag_posi
EZH2_CBX2.methy %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  tidyr::nest(-tag,-group,-Gene_Symbol) %>%
  dplyr::group_by(tag,group) %>%
  dplyr::mutate(median=purrr::map(data,fn_median)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(EZH2_CBX2.tag_posi,by="tag") %>%
  dplyr::ungroup()-> EZH2_CBX2.median.methy

# calculation -------------------------------------------------------------

# EZH2 and CBX2 -----------------------------------------------------------
EZH2_CBX2.methy %>%
  tidyr::drop_na(methy) %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  dplyr::group_by(tag) %>%
  dplyr::do(
    broom::tidy(
      t.test(methy ~ group, data = .)
    )
  ) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.01,"**",sig)) %>%
  dplyr::select(tag,p.value,sig,alternative)  -> EZH2_CBX2.ttest

EZH2_CBX2.methy %>%
  dplyr::inner_join(EZH2_CBX2.tag_posi,by="tag") %>%
  dplyr::mutate(methy=as.numeric(methy)) %>%
  dplyr::inner_join(EZH2_CBX2.ttest,by="tag") %>%
  dplyr::mutate(lab.y=1.05) %>%
  ggplot(aes(x=Genomic_Direction,y=methy,color=group)) +
  geom_point(position = "jitter",size=1) +
  geom_boxplot() +
  geom_text(aes(x=Genomic_Direction,y=lab.y,label=sig),color="black") +
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
  ) -> EZH2_CBX2.box;EZH2_CBX2.box
ggsave(file.path(out_path_fig,"Figure3","Figure3B.Methy_histone.pdf"),EZH2_CBX2.box,device = "pdf",width = 10,height = 6)

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
  dplyr::group_by(tag) %>%
  dplyr::do(
    broom::tidy(
      t.test(methy ~ group, data = .)
    )
  ) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(p.value<=0.01,"**",sig)) %>%
  dplyr::mutate(diff = -estimate) %>%
  dplyr::select(tag,p.value,diff) -> genelist_pro.ttest

genelist_pro.ttest %>%
  dplyr::inner_join(genelist_pro.tag_posi,by="tag") %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::mutate(diff.m=mean(diff)) %>%
  dplyr::select(Gene_Symbol,diff.m) %>%
  unique()

