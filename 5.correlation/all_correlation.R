
# .libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(magrittr)
library(corrplot)
# data path ---------------------------------------------------------------

miRNA_list_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"


gene_list_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
# output data -------------------------------------------------------------
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# laod data ---------------------------------------------------------------

miRNA_list <- readr::read_tsv(file.path(miRNA_list_path,"NOISeq_DE_mirna_FC2_cpm30.mirnaid.up.id"))
mirna_fc1.5_cpm5000 <- readr::read_tsv(file.path("H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/no_FC/NOISeq_DE_mirna_noFC_cpm30.mirnaid")) %>%
  dplyr::filter(abs(log2FC)>=0.585) %>%
  dplyr::filter(case_mean>5000|con_mean>5000) %>%
  dplyr::select(mirna)
miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
genelist_TF <- readr::read_tsv(file.path(gene_list_path,"all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(Class=="Down")
# genelist_pro <- readr::read_tsv(file.path(gene_list_path,"DOWN1.5_exp30_pro_EHZ2_CBX2_common_targets.DE_info"))

# data manage -------------------------------------------------------------
miRNA_exp %>%
  dplyr::filter(name %in% c(miRNA_list$mirna_id,mirna_fc1.5_cpm5000$mirna)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> miRNA_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(genelist_pro$gene_id.x)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> gene_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c("EZH2","CBX2")) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> EZH2_CBX2_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c("E2F1","SOX9","E2F3","FOXP3","SOX4")) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> EZH2_CBX2_upstreamTF_exp.gather

# CBX2_upstreamMIR<- c("hsa-miR-331-3p","hsa-let-7f-5p","hsa-miR-93-5p","hsa-miR-629-5p","hsa-miR-338-3p","hsa-miR-1-3p","hsa-miR-30a-5p",
#   "hsa-miR-193b-3p","hsa-miR-128-3p","hsa-miR-183-5p","hsa-miR-1287-5p","hsa-let-7c-5p","hsa-miR-195-5p","hsa-miR-424-5p",
#   "hsa-miR-2355-5p")
# EZH2_upstreamMIR<- c("hsa-miR-200a-3p","hsa-miR-139-5p","hsa-miR-93-5p","hsa-miR-193b-3p","hsa-miR-128-3p","hsa-let-7c-5p",
#                      "hsa-miR-195-5p","hsa-miR-429","hsa-miR-92b-3p")
EZH2_CBX2_common_miR <- c("hsa-let-7a-5p","hsa-let-7c-5p","hsa-miR-30d-5p","hsa-miR-101-3p","hsa-miR-195-5p")

# miRNA_exp %>%
#   dplyr::filter(name %in% CBX2_upstreamMIR) %>%
#   dplyr::select(-cancer_types,-gene) %>%
#   tidyr::gather(-name,key="sample",value="mirna_exp") %>%
#   dplyr::mutate(sample=substr(sample,9,16)) %>%
#   dplyr::as_tibble() %>%
#   tidyr::nest(-name) -> CBX2_upstreamMIR_exp.gather
# miRNA_exp %>%
#   dplyr::filter(name %in% EZH2_upstreamMIR) %>%
#   dplyr::select(-cancer_types,-gene) %>%
#   tidyr::gather(-name,key="sample",value="mirna_exp") %>%
#   dplyr::mutate(sample=substr(sample,9,16)) %>%
#   dplyr::as_tibble() %>%
#   tidyr::nest(-name) -> EZH2_upstreamMIR_exp.gather
miRNA_exp %>%
  dplyr::filter(name %in% EZH2_CBX2_common_miR) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> EZH2_CBX2_upstreamMIR_exp.gather
# calculation -------------------------------------------------------------

# miRNA and EZH2&CBX2 targets ----
#funciton ----
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

miRNA_exp.gather %>%
  # dplyr::filter(name %in% c("hsa-miR-21-5p")) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=gene_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> mirna_gene_spm_cor
mirna_gene_spm_cor %>%
  readr::write_tsv(file.path(out_path_sup,"correlation-table","mirna_EZH2-CBX2-targets_spm_cor.tsv"))
mirna_gene_spm_cor %>%
  dplyr::select(name,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> gene_Cor_matrix
rownames(gene_Cor_matrix) <- gene_Cor_matrix$name
gene_Cor_matrix <- gene_Cor_matrix[,-1]

mirna_gene_spm_cor %>%
  dplyr::mutate(Sig=ifelse(Cor < 0 & p.value<=0.05,"*","")) %>%
  dplyr::mutate(Sig=ifelse(Cor < -0.4 & p.value<=0.05,"**",Sig)) %>%
  dplyr::select(name,symbol,Sig) %>%
  tidyr::spread(symbol,Sig) %>%
  as.data.frame()-> gene_sig_matrix
rownames(gene_sig_matrix) <- gene_sig_matrix$name
gene_sig_matrix <- gene_sig_matrix[,-1]
library(pheatmap)

pdf(file.path(out_path_sup,"Supplementary-Fig4.EZH2-CBX2_upFC2-miRNA_spearman.pdf"),height = 8,width = 15)
pheatmap(gene_Cor_matrix,color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(30),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,cutree_cols = 2,
         treeheight_row=0,treeheight_col = 0,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10,legend = TRUE,
         display_numbers = gene_sig_matrix)
dev.off()

# EZH2&CBX2 with their targets ----

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
EZH2_CBX2_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=gene_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> EZH2_CBX2_targets_spm_cor
EZH2_CBX2_targets_spm_cor %>%
  readr::write_tsv(file.path(out_path_sup,"correlation-table","EZH2_CBX2_targets_spm_cor.tsv"))

EZH2_CBX2_targets_spm_cor %>%
  dplyr::select(symbol1,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> EZH2_CBX2_targets_spm_cor_matrix
rownames(EZH2_CBX2_targets_spm_cor_matrix) <- EZH2_CBX2_targets_spm_cor_matrix$symbol1
EZH2_CBX2_targets_spm_cor_matrix <- EZH2_CBX2_targets_spm_cor_matrix[,-1]

EZH2_CBX2_targets_spm_cor %>%
  dplyr::mutate(Sig=ifelse(Cor < 0 & p.value<=0.05,"*","")) %>%
  dplyr::mutate(Sig=ifelse(Cor < -0.4 & p.value<=0.05,"**",Sig)) %>%
  dplyr::select(symbol1,symbol,Sig) %>%
  tidyr::spread(symbol,Sig) %>%
  as.data.frame() -> EZH2_CBX2_targets_sig_matrix
rownames(EZH2_CBX2_targets_sig_matrix) <- EZH2_CBX2_targets_sig_matrix$symbol1
EZH2_CBX2_targets_sig_matrix <- EZH2_CBX2_targets_sig_matrix[,-1]

pdf(file.path(out_path_sup,"Supplementary-Fig5.EZH2-CBX2_targets_spearman.pdf"),height = 10,width = 4)
pheatmap(head(EZH2_CBX2_targets_spm_cor_matrix,44),color = colorRampPalette(c("deepskyblue3","white"))(30),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=0,treeheight_col = 0,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 8,cellheight = 8,legend = TRUE,
         display_numbers = head(EZH2_CBX2_targets_sig_matrix,44))
pheatmap(EZH2_CBX2_targets_spm_cor_matrix[45:88,],color = colorRampPalette(c("deepskyblue3","white"))(30),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=0,treeheight_col = 0,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 8,cellheight = 8,legend = TRUE,
         display_numbers = EZH2_CBX2_targets_sig_matrix[45:88,])
dev.off()

# EZH2&CBX2 with their upstream TF ----
EZH2_CBX2_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=EZH2_CBX2_upstreamTF_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> EZH2_CBX2_upstreamTF_spm_cor
EZH2_CBX2_upstreamTF_spm_cor %>%
  readr::write_tsv(file.path(out_path_sup,"correlation-table","EZH2_CBX2_upstreamTF_spm_cor.tsv"))
#pic drawing
EZH2_CBX2_upstreamTF_spm_cor %>%
  dplyr::select(symbol1,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> EZH2_CBX2_upstreamTF_spm_cor_matrix
rownames(EZH2_CBX2_upstreamTF_spm_cor_matrix) <- EZH2_CBX2_upstreamTF_spm_cor_matrix$symbol1
EZH2_CBX2_upstreamTF_spm_cor_matrix <- EZH2_CBX2_upstreamTF_spm_cor_matrix[,-1]

EZH2_CBX2_upstreamTF_spm_cor %>%
  dplyr::mutate(Sig=ifelse(Cor >= 0.3 & p.value<=0.05,"*","")) %>%
  dplyr::select(symbol1,symbol,Sig) %>%
  tidyr::spread(symbol,Sig) %>%
  as.data.frame() -> EZH2_CBX2_upstreamTF_sig_matrix
rownames(EZH2_CBX2_upstreamTF_sig_matrix) <- EZH2_CBX2_upstreamTF_sig_matrix$symbol1
EZH2_CBX2_upstreamTF_sig_matrix <- EZH2_CBX2_upstreamTF_sig_matrix[,-1]

pdf(file.path(out_path_fig,"Figure2","Figure2C.EZH2-CBX2_upstream_TFs_spearman.pdf"),height = 4,width = 4)
pheatmap(EZH2_CBX2_upstreamTF_spm_cor_matrix,color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(20),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=0,treeheight_col = 0,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10, legend = TRUE,
         display_numbers = EZH2_CBX2_upstreamTF_sig_matrix)
dev.off()

# EZH2 and CBX2 with common upstream miRs ------
EZH2_CBX2_upstreamMIR_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp= EZH2_CBX2_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> CBX2_EZH2_upstreamMIR_spm_cor

CBX2_EZH2_upstreamMIR_spm_cor %>%
  readr::write_tsv(file.path(out_path_sup,"correlation-table","CBX2_EZH2_upstreamMIR_spm_cor.tsv"))
#pic drawing
CBX2_EZH2_upstreamMIR_spm_cor %>%
  dplyr::select(name,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> CBX2_EZH2_upstreamMIR_spm_cor_matrix
rownames(CBX2_EZH2_upstreamMIR_spm_cor_matrix) <- CBX2_EZH2_upstreamMIR_spm_cor_matrix$name
CBX2_EZH2_upstreamMIR_spm_cor_matrix <- CBX2_EZH2_upstreamMIR_spm_cor_matrix[,-1]

CBX2_EZH2_upstreamMIR_spm_cor %>%
  dplyr::mutate(Sig=ifelse(Cor <= -0.3 & p.value<=0.05,"*","")) %>%
  dplyr::select(name,symbol,Sig) %>%
  tidyr::spread(symbol,Sig) %>%
  as.data.frame() -> CBX2_EZH2_upstreamMIR_sig_matrix
rownames(CBX2_EZH2_upstreamMIR_sig_matrix) <- CBX2_EZH2_upstreamMIR_sig_matrix$name
CBX2_EZH2_upstreamMIR_sig_matrix <- CBX2_EZH2_upstreamMIR_sig_matrix[,-1]

pdf(file.path(out_path_fig,"Figure2","Figure2C.EZH2-CBX2_upstream_miRs_spearman.pdf"),height = 4,width = 4)
pheatmap(rbind(CBX2_EZH2_upstreamMIR_spm_cor_matrix,EZH2_CBX2_upstreamTF_spm_cor_matrix),color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(20),
         kmeans_k = NA,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=0,treeheight_col = 0,breaks = seq(-1,1,0.1),
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10, legend = TRUE,
         display_numbers = rbind(CBX2_EZH2_upstreamMIR_sig_matrix,EZH2_CBX2_upstreamTF_sig_matrix))
dev.off()

# CBX2 with upstram miRs ----
EZH2_CBX2_exp.gather %>%
  dplyr::filter(symbol=="CBX2") -> CBX2_exp.gather
CBX2_upstreamMIR_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=CBX2_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> CBX2_upstreamMIR_spm_cor

CBX2_upstreamMIR_spm_cor %>%
  dplyr::select(name,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> CBX2_upstreamMIR_spm_cor_matrix
cnx2_mir_name <- CBX2_upstreamMIR_spm_cor_matrix$name
CBX2_upstreamMIR_spm_cor_matrix <- as.matrix(CBX2_upstreamMIR_spm_cor_matrix[,-1],nrow=1)
rownames(CBX2_upstreamMIR_spm_cor_matrix) <- cnx2_mir_name
colnames(CBX2_upstreamMIR_spm_cor_matrix) <- "CBX2"

CBX2_upstreamMIR_spm_cor %>%
  dplyr::mutate(Sig=ifelse(Cor < 0 & p.value<=0.05,"*","")) %>%
  dplyr::mutate(Sig=ifelse(Cor < -0.4 & p.value<=0.05,"**",Sig)) %>%
  dplyr::select(name,symbol,Sig) %>%
  tidyr::spread(symbol,Sig) %>%
  as.data.frame() -> CBX2_upstreamMIR_sig_matrix
cnx2_mir_name_1 <- CBX2_upstreamMIR_sig_matrix$name
CBX2_upstreamMIR_sig_matrix <- as.matrix(CBX2_upstreamMIR_sig_matrix[,-1],nrow=1)
rownames(CBX2_upstreamMIR_sig_matrix) <- cnx2_mir_name
colnames(CBX2_upstreamMIR_sig_matrix) <- "CBX2"

pdf(file.path(out_path_fig,"Figure2","Figure2C.CBX2_upstream_miRs_spearman.pdf"),height = 4,width = 4)
pheatmap(CBX2_upstreamMIR_spm_cor_matrix,color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(20),
         kmeans_k = NA,cluster_cols = FALSE,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=30,treeheight_col = 30,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10, legend = TRUE,
         display_numbers = CBX2_upstreamMIR_sig_matrix)
dev.off()

# EZH2 with upstram miRs ----
EZH2_CBX2_exp.gather %>%
  dplyr::filter(symbol=="EZH2") -> EZH2_exp.gather
EZH2_upstreamMIR_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=EZH2_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> EZH2_upstreamMIR_spm_cor

EZH2_upstreamMIR_spm_cor %>%
  dplyr::select(name,symbol,Cor) %>%
  tidyr::spread(symbol,Cor) %>%
  as.data.frame() -> EZH2_upstreamMIR_spm_cor.matrix
ezh2_mir_name <- EZH2_upstreamMIR_spm_cor.matrix$name
EZH2_upstreamMIR_spm_cor.matrix <- as.matrix(EZH2_upstreamMIR_spm_cor.matrix[,-1],ncol=1)
dimnames(EZH2_upstreamMIR_spm_cor.matrix)=list(ezh2_mir_name,c("EZH2"))

EZH2_upstreamMIR_spm_cor %>%
  dplyr::mutate(sig=ifelse(Cor < 0 & p.value<=0.05,"*","")) %>%
  dplyr::mutate(sig=ifelse(Cor < -0.4 & p.value<=0.05,"**",sig)) %>%
  dplyr::select(name,symbol,sig) %>%
  tidyr::spread(symbol,sig) %>%
  as.data.frame() -> EZH2_upstreamMIR_sig.matrix
EZH2_upstreamMIR_sig.matrix <- as.matrix(EZH2_upstreamMIR_sig.matrix[,-1],ncol=1)
dimnames(EZH2_upstreamMIR_sig.matrix)=list(ezh2_mir_name,c("EZH2"))

pdf(file.path(out_path_fig,"Figure2","Figure2C.EZH2_upstream_miRs_spearman.pdf"),height = 4,width = 4)
pheatmap(EZH2_upstreamMIR_spm_cor.matrix,color = colorRampPalette(c("deepskyblue3","white","lightcoral"))(20),
         kmeans_k = NA,cluster_cols = FALSE,border_color = "grey60",cutree_rows = 2,#cutree_cols = 2,
         treeheight_row=30,treeheight_col = 30,
         fontsize =15,fontsize_row = 10,fontsize_col = 10, cellwidth = 10,cellheight = 10, legend = TRUE,
         display_numbers = EZH2_upstreamMIR_sig.matrix)
dev.off()
# draw pic ----------------------------------------------------------------

# up miR and EZH2 & CBX2 targtes ----


mirna_gene_spm_cor %>%
  readr::write_tsv(file.path(out_path,"Supplementary_Table2.up_FC2_mirna-EZH2_CBX2_targets-spearman.xls"))

################ cell cycle pathway regulation pairs -----

cellcycle_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-cellcycle-relatedgenes"
cellcycle_genelist <- readr::read_tsv(file.path(cellcycle_path,"attribute.hallmark-added.txt"))
cellcycle_TF <- cellcycle_genelist %>% dplyr::filter(gene_type==1)
cellcycle_miR <- cellcycle_genelist %>% dplyr::filter(gene_type==2)
cellcycle_gene <- cellcycle_genelist %>% dplyr::filter(gene_type==3)

miRNA_exp %>%
  dplyr::filter(name %in% c(cellcycle_miR$SYMBOL)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> cellcycle_miRNA_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(cellcycle_TF$SYMBOL)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> cellcycle_TF_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(cellcycle_gene$SYMBOL,cellcycle_TF$SYMBOL)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> cellcycle_gene_exp.gather

cellcycle_miRNA_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=cellcycle_gene_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="name","gene"="symbol") -> cellcycle_miR_gene_spm_cor

cellcycle_TF_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=cellcycle_gene_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="symbol","gene"="symbol1") -> cellcycle_TF_gene_spm_cor


cellcyle_net <- readr::read_tsv(file.path(cellcycle_path,"network.txt")) %>%
  dplyr::rename("regulator"="From") %>%
  dplyr::filter(regulate_type==1 | regulate_type==3)
cellcycle_miR_gene_spm_cor %>%
  rbind(cellcycle_TF_gene_spm_cor) %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(cellcyle_net,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene==To) %>%
  dplyr::select(regulator,gene,Cor,p.value)-> cellcyle_all_cor

# draw pic ----
library(ggplot2)
library(guide)
CPCOLS <- c("red", "white", "blue")
cellcyle_all_cor %>%
  dplyr::filter(p.value<=0.05 & abs(Cor)>0.3) %>%
  dplyr::mutate(`-log10(P)`=-log10(p.value)) %>%
  dplyr::mutate(`-log10(P)`=ifelse(`-log10(P)`=="Inf" |`-log10(P)`>5,5,`-log10(P)`)) -> ready_draw
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

ready_draw %>%
  ggplot(aes(x=gene,y=regulator,color=Cor)) +
  geom_tile(fill=c("#EDEDED"),colour = "grey") +
  geom_point(size=3) +
  guides(color=guide_colorbar(title.position="left")) +
  scale_x_discrete(limits = x_rank$gene) +
  scale_y_discrete(limits = y_rank$regulator) +
  ylab("Regulators") +
  xlab("Cell cycle genes") +
  scale_color_gradient2(
    name = "Spearman r", # "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1],
    breaks=c(-0.4,0,0.4,0.8)
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
ggsave(file.path(out_path_fig,"Figure1","Figue S2C.cell_cycle.correlation.pdf"),device = "pdf",width = 10,height = 4)
ggsave(file.path(out_path_fig,"Figure1","Figue S2C.cell_cycle.correlation.tiff"),device = "tiff",width = 10,height = 4)


########### PPAR pathway regulation pairs -----

ppar_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea/FFL/LUAD-noFC-prob0.9-kegg-gsea-ppar-relatedgenes"
ppar_genelist <- readr::read_tsv(file.path(ppar_path,"attribute.hallmark-added.txt"))
ppar_TF <- ppar_genelist %>% dplyr::filter(gene_type==1)
ppar_miR <- ppar_genelist %>% dplyr::filter(gene_type==2)
ppar_gene <- ppar_genelist %>% dplyr::filter(gene_type==3)

miRNA_exp %>%
  dplyr::filter(name %in% c(ppar_miR$SYMBOL)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="mirna_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-name) -> ppar_miRNA_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(ppar_TF$SYMBOL)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> ppar_TF_exp.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(ppar_gene$SYMBOL,ppar_TF$SYMBOL)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::filter(substr(sample,6,6)=="0") %>%
  dplyr::as_tibble() %>%
  tidyr::nest(-symbol) -> ppar_gene_exp.gather

ppar_miRNA_exp.gather %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=ppar_gene_exp.gather,fn_get_spm_a)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="name","gene"="symbol") -> ppar_miR_gene_spm_cor

ppar_TF_exp.gather %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(spm=purrr::map(data,g_exp=ppar_gene_exp.gather,fn_get_spm_a_g)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::rename("regulator"="symbol","gene"="symbol1") -> ppar_TF_gene_spm_cor


ppar_net <- readr::read_tsv(file.path(ppar_path,"network.txt")) %>%
  dplyr::rename("regulator"="From") %>%
  dplyr::filter(regulate_type==1 | regulate_type==3)
ppar_miR_gene_spm_cor %>%
  rbind(ppar_TF_gene_spm_cor) %>%
  dplyr::select(regulator,gene,Cor,p.value) %>%
  dplyr::right_join(ppar_net,by="regulator") %>%
  tidyr::drop_na() %>%
  dplyr::filter(gene==To) %>%
  dplyr::select(regulator,gene,Cor,p.value)-> ppar_all_cor

# draw pic ----
library(ggplot2)

CPCOLS <- c("red", "white", "blue")
ppar_all_cor %>%
  dplyr::filter(p.value<=0.05 & abs(Cor)>=0.3) %>%
  dplyr::mutate(`-log10(P)`=-log10(p.value)) %>%
  dplyr::mutate(`-log10(P)`=ifelse(`-log10(P)`=="Inf" |`-log10(P)`>5,5,`-log10(P)`)) -> ready_draw
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

ready_draw %>%
  ggplot(aes(x=gene,y=regulator,color=Cor)) +
  geom_tile(fill=c("#EDEDED"),colour = "grey") +
  geom_point(size=5) +
  guides(color=guide_colorbar(title.position="left")) +
  scale_x_discrete(limits = x_rank$gene) +
  scale_y_discrete(limits = y_rank$regulator) +
  ylab("Regulators") +
  xlab("PPAR signaling pathway genes") +
  scale_color_gradient2(
    name = "Spearman r", # "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1]
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
    legend.position = c(0.4,0.3),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20)
  ) -> p;p
ggsave(file.path(out_path_fig,"Figure1","Figue S2D.ppar.correlation.pdf"),device = "pdf",width = 3,height = 3)
ggsave(file.path(out_path_fig,"Figure1","Figue S2D.ppar.correlation.tiff"),device = "tiff",width = 3,height = 3)
