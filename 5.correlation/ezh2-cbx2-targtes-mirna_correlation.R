
.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(magrittr)
# data path ---------------------------------------------------------------

miRNA_list_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
miRNA_exp_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/肺癌基因表达量数据/rsem标准化表达量"
gene_list_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"

# laod data ---------------------------------------------------------------

miRNA_list <- readr::read_tsv(file.path(miRNA_list_path,"NOISeq_DE_mirna_FC2_cpm30.mirnaid.up.id"))
miRNA_exp <- read.table(file.path(miRNA_exp_path,"miRNA_isoform.expression.mirnaid.noMIMA"),sep="\t",header = T)
pro_exp <- read.table(file.path(miRNA_exp_path,"Pro_rsem_exp.xls"),sep="\t",header = T)
TF_exp <- read.table(file.path(miRNA_exp_path,"TF_rsem_exp.xls"),sep="\t",header = T)
genelist_pro <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
genelist_TF <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))

# data manage -------------------------------------------------------------
miRNA_exp %>%
  dplyr::filter(mirna_id %in% miRNA_list$mirna_id) %>%
  tidyr::gather(-mirna_id,key="sample",value="mirna_exp") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-mirna_id) -> miRNA_exp.gather

pro_exp %>%
  rbind(TF_exp) %>%
  dplyr::filter(gene_id %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) %>%
  tidyr::gather(-gene_id,key="sample",value="gene_exp") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-gene_id) -> gene_exp.gather

fn_get_spm_a <- function(m_data){
  # print(m_data)
  gene_exp.gather %>%
    dplyr::group_by(gene_id) %>%
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
miRNA_exp.gather$mirna_id <- as.character(miRNA_exp.gather$mirna_id)
gene_exp.gather$gene_id <- as.character(gene_exp.gather$gene_id)
miRNA_exp.gather %>%
  # dplyr::filter(mirna_id %in% c("hsa-miR-21-5p")) %>%
  dplyr::group_by(mirna_id) %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> mirna_gene_spm_cor

mirna_gene_spm_cor %>%
  dplyr::select(mirna_id,gene_id,Cor) %>%
  tidyr::spread(gene_id,Cor) %>%
  as.data.frame() -> gene_Cor_matrix
rownames(gene_Cor_matrix) <- gene_Cor_matrix$mirna_id
gene_Cor_matrix <- gene_Cor_matrix[,-1]

mirna_gene_spm_cor %>%
  dplyr::select(mirna_id,gene_id,p.value) %>%
  tidyr::spread(gene_id,p.value) -> gene_pvalue_matrix
library(pheatmap)

pheatmap(gene_Cor_matrix,color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(30),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,cutree_cols = 2,
         treeheight_row=30,treeheight_col = 30,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10,legend = TRUE,
         display_numbers = matrix(ifelse(gene_Cor_matrix <= (-0.4), "*", ""),nrow(gene_Cor_matrix)))

         
# output data -------------------------------------------------------------
out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
mirna_gene_spm_cor %>%
  readr::write_tsv(file.path(out_path,"Supplementary_Table2.up_FC2_mirna-EZH2_CBX2_targets-spearman.xls"))
