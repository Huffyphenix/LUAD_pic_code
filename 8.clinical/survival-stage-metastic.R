
# data path ---------------------------------------------------------------
target_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
exp_path <- "H:/data/TCGA/TCGA_data"
clinical_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/生存分析/data/LUAD"
survival_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure5"
# load data ---------------------------------------------------------------

targets_tf <- readr::read_tsv(file.path(target_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.TF_LUAD-FC2.down"))
target_pro <- readr::read_tsv(file.path(target_path,"CBX2-pva3_H3K27me3-pva3_OVERLAP-100bp_GRCh38-hg38_TSS-5kb.gene_symbol.protein_coding_LUAD-FC2_down"))
miRNA_exp <- readr::read_rds(file.path(exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()
clinical <- readr::read_tsv(file.path(clinical_path,"clinical_info_multi_survival.txt"))
EZH2upstream_miRNA_list <- c("hsa-miR-101-3p","hsa-miR-30d-5p","hsa-let-7c-5p")
targetsupstream_miRNA_list <- c("hsa-miR-210-3p","hsa-miR-183-5p","hsa-miR-151a-5p","hsa-miR-1307-3p","hsa-miR-93-5p",
                                "hsa-miR-2355-5p","hsa-miR-141-5p","hsa-miR-141-92b-3p","hsa-miR-141-3p")
keys <- c("EZH2","CBX2","E2F1","SOX4")


# filter ------------------------------------------------------------------

miRNA_exp %>%
  dplyr::filter(name %in% c(EZH2upstream_miRNA_list,targetsupstream_miRNA_list)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==0) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> miRNA_exp.gather
miRNA_exp %>%
  dplyr::filter(name %in% c(EZH2upstream_miRNA_list,targetsupstream_miRNA_list)) %>%
  dplyr::select(-cancer_types,-gene) %>%
  tidyr::gather(-name,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==1) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> miRNA_exp_normal.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(targets_tf$gene_id,target_pro$gene_id,keys)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==0) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> mRNA_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(targets_tf$gene_id,target_pro$gene_id,keys)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==1) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> mRNA_exp_normal.gather

clinical %>%
  dplyr::mutate(status=ifelse(vital_status=="dead",as.integer(0),as.integer(1))) %>%
  dplyr::mutate(OS=ifelse(is.na(days_to_death),days_to_last_followup,days_to_death)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage i","Stage I","stage III")) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ia","Stage I",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ib","Stage I",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ii","Stage II",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage iia","Stage II",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage iv","Stage IV",stage)) %>%
  dplyr::select(patient_id,status,OS,stage,pathologic_m) %>%
  dplyr::rename("sample"="patient_id")-> clinical_data

                

# Survival ----------------------------------------------------------------

#### miRNA ----
# function----
fn_survival <- function(.gene,.data,.up,.low){
  # .up=50
  # .low=50
  # .gene = "xxx"
  .data %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> .data
  # print(.d)
  .d_diff <- survival::survdiff(survival::Surv(OS, status) ~ group, data = .data)
  # print(.d_diff)
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
  # print(kmp)
  coxp <-  broom::tidy(survival::coxph(survival::Surv(OS, status) ~ exp, data = .data, na.action = na.exclude))
  # print(coxp)
  
  .fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = .data, na.action = na.exclude) 
  # print(.fit)
  purrr::pwalk(.gene,kmp,.data,fn_sur_draw)
  data.frame(KMP=kmp,Coxp=coxp) -> .out
  return(.out)
}
fn_sur_draw <- function(.gene,kmp,.data){
  # miRNA_clinical %>%
  #   dplyr::filter(name==name) %>%
  #   tidyr::unnest() %>%
  #   dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
  #   dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) ->  .d
  
  fig_name <- paste(.gene,kmp,"pdf", sep = ".")
  .fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = .data, na.action = na.exclude) 
  survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                          surv.median.line = "hv",
                          title = paste(paste(name, sep = "-"), "Logrank P =", signif(KMP, 2)),
                          xlab = "Survival in days",
                          ylab = 'Probability of survival',
                          legend.title = "Expression group:",
                          legend= c(0.8,0.8),
                          risk.table = TRUE,
                          tables.height = 0.2,
                          palette = c("#E7B800", "#2E9FDF"),
                          ggtheme = theme_bw()
    )
  ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path), width = 6, height = 6)
}
  
miRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-name) -> miRNA_clinical
miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(survival=purrr::map2(name,data,.up=50,.low=50,fn_survival)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> miRNA_survival_p
miRNA_survival_p %>%
  dplyr::group_by(name) %>%
  dplyr::select(name,KMP) %>%
  purrr::pwalk(miRNA_clinical=miRNA_clinical,,.up=50,.low=50,.f=fn_sur_draw)
