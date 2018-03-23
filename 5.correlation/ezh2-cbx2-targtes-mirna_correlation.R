
.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(magrittr)
# data path ---------------------------------------------------------------

miRNA_list_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
miRNA_exp_path <- "F:/data/TCGA/TCGA_data"
gene_list_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"

# laod data ---------------------------------------------------------------

miRNA_list <- readr::read_tsv(file.path(miRNA_list_path,"NOISeq_DE_mirna_FC2_cpm30.mirnaid.up.id"))
miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
genelist_pro <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
genelist_TF <- readr::read_tsv(file.path(gene_list_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))

# data manage -------------------------------------------------------------
miRNA_exp %>%
  dplyr::filter(name %in% miRNA_list$mirna_id) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(9,16,sample)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> miRNA_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(genelist_pro$gene_id,genelist_TF$gene_id)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(9,16,sample)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> gene_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c("EZH2","CBX2")) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(9,16,sample)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> EZH2_CBX2_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c("E2F1","SOX9","EGR2","EGR1","NEGR1","E2F3","KLF6","FOXP3","SOX4")) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(9,16,sample)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> EZH2_CBX2_upstreamTF_exp.gather

CBX2_upstreamMIR<- c("hsa-miR-331-3p","hsa-let-7f-5p","hsa-miR-93-5p","hsa-miR-629-5p","hsa-miR-338-3p","hsa-miR-1-3p","hsa-miR-30a-5p",
  "hsa-miR-193b-3p","hsa-miR-128-3p","hsa-miR-183-5p","hsa-miR-1287-5p","hsa-let-7c-5p","hsa-miR-195-5p","hsa-miR-424-5p",
  "hsa-miR-2355-5p")
EZH2_upstreamMIR<- c("hsa-miR-200a-3p","hsa-miR-139-5p","hsa-miR-93-5p","hsa-miR-193b-3p","hsa-miR-128-3p","hsa-let-7c-5p",
                     "hsa-miR-195-5p","hsa-miR-429","hsa-miR-92b-3p")
miRNA_exp %>%
  dplyr::filter(name %in% CBX2_upstreamMIR) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> CBX2_upstreamMIR_exp.gather
miRNA_exp %>%
  dplyr::filter(name %in% EZH2_upstreamMIR) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> EZH2_upstreamMIR_exp.gather
# calculation -------------------------------------------------------------
#funciton ----
fn_get_spm_a <- function(m_data){
  # print(m_data)
  gene_exp.gather %>%
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
miRNA_exp.gather$mirna_id <- as.character(miRNA_exp.gather$mirna_id)
gene_exp.gather$gene_id <- as.character(gene_exp.gather$gene_id)

# miRNA and EZH2&CBX2 targets
miRNA_exp.gather %>%
  dplyr::filter(name %in% c("hsa-miR-21-5p")) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> mirna_gene_spm_cor

# EZH2&CBX2 with their targets
gene_exp.gather %>%
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
