.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(magrittr,ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
# data path ---------------------------------------------------------------

# target_path <- "S:??????/?ҵļ?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
target_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
exp_path <- "H:/data/TCGA/TCGA_data"
clinical_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/生存分析/data/LUAD"
clinical_path_1 <- "F:/我的坚果云/ENCODE-TCGA-LUAD/survival"
survival_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure5"
data_path<- "H:/data"
chip_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/"

# survival_path <- "S:??????/?ҵļ?????/ENCODE-TCGA-LUAD/Figure/Figure5"
# load data ---------------------------------------------------------------
Animal_TF <-  readr::read_tsv(file.path(data_path,"AnimalTFDB","Homo_sapiens_transcription_factors_gene_list.txt"))
enzyme_lsit <- readr::read_tsv(file.path(chip_path,"enzyme_list.symbol.xls")) %>%
  .$Symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Animal_TF %>%
  dplyr::filter(! Entrez_ID %in% enzyme_lsit$ENTREZID) -> Animal_TF
Animal_TF$Entrez_ID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db) -> Animal_TF.symbol

targets <- readr::read_tsv(file.path(target_path,"all_EHZ2_CBX2_common_targets.DE_info"))
targets %>%
  dplyr::filter(prob>0.99 & log2FC <= -0.585) %>%
  dplyr::filter(case_mean>=30) %>%
  dplyr::filter(gene_id.x %in% Animal_TF.symbol$SYMBOL) -> targets_down_TF
targets_down_TF %>%
  readr::write_tsv(file.path(target_path,"DOWN1.5_allexp30_TF_EHZ2_CBX2_common_targets.DE_info"))
targets %>%
  dplyr::filter(prob>0.99 & log2FC <= -0.585) %>%
  dplyr::filter(case_mean>=30) %>%
  dplyr::filter(! gene_id.x %in% Animal_TF.symbol$SYMBOL) -> targets_down_pro
targets_down_pro %>%
  readr::write_tsv(file.path(target_path,"DOWN1.5_allexp30_pro_EHZ2_CBX2_common_targets.DE_info"))

miRNA_exp <- readr::read_rds(file.path(exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()
clinical <- readr::read_tsv(file.path(clinical_path,"clinical_info_multi_survival.txt"))
survival <- readr::read_tsv(file.path(clinical_path_1,"cancer_cell_survival_time_table.txt"))
survival %>%
  dplyr::mutate(sample=substr(bcr_patient_barcode,9,12)) %>%
  dplyr::select(sample,PFI.1,PFI.time.1) %>%
  dplyr::mutate(PFI.time.1=as.numeric(PFI.time.1)) -> PFI_survival_time
  
EZH2upstream_miRNA_list <- c("hsa-miR-101-3p","hsa-miR-30d-5p","hsa-let-7c-5p")
targetsupstream_miRNA_list <- c("hsa-miR-210-3p","hsa-miR-183-5p","hsa-miR-151a-5p","hsa-miR-1307-3p","hsa-miR-93-5p",
                                "hsa-miR-2355-5p","hsa-miR-141-5p","hsa-miR-92b-3p","hsa-miR-141-3p","hsa-miR-192-5p",
                                "hsa-miR-194-5p","hsa-let-7i-3p","hsa-miR-106b-3p","hsa-miR-424-5p")
keys <- c("EZH2","CBX2","E2F1","SOX4")


# filter ------------------------------------------------------------------
library(magrittr)
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
  dplyr::as_tibble() %>%
  dplyr::mutate(stage="Normal(TA)") -> miRNA_exp_normal.gather

gene_exp %>%
  dplyr::filter(symbol %in% c(targets_down_TF$gene_id.x,targets_down_pro$gene_id.x,keys)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==0) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> mRNA_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(targets_down_TF$gene_id.x,targets_down_pro$gene_id.x,keys)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==1) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(stage="Normal(TA)") -> mRNA_exp_normal.gather

clinical %>%
  dplyr::mutate(status=ifelse(vital_status=="dead",as.integer(1),as.integer(0))) %>%
  dplyr::mutate(OS=ifelse(is.na(days_to_death),days_to_last_followup,days_to_death)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage i","stage I","stage III")) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ia","stage I",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ib","stage I",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage ii","stage II",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage iia","stage II",stage)) %>%
  dplyr::mutate(stage=ifelse(pathologic_stage=="stage iv","stage IV",stage)) %>%
  dplyr::mutate(metastasis=ifelse(pathologic_m=="m0","m0","m1")) %>%
  dplyr::mutate(metastasis=ifelse(pathologic_m=="mx",NA,metastasis)) %>%
  dplyr::select(patient_id,status,OS,stage,metastasis) %>%
  dplyr::rename("sample"="patient_id")-> clinical_data

clinical_data %>%
  dplyr::full_join(PFI_survival_time,by="sample") -> clinical_data

# Survival ----------------------------------------------------------------

#### miRNA ----
# function----
fn_survival <- function(.gene,.data,.up,.low){
  # .up=50
  # .low=50
  # .gene = "xxx"
  print(.gene)
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
  # fn_sur_draw(.gene,kmp,.fit)
  data.frame(KMP=kmp,Coxp=coxp) -> .out
  return(.out)
}
fn_survival_PFI <- function(.gene,.data,.up,.low){
  # .up=50
  # .low=50
  # .gene = "xxx"
  .data %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> .data
  # print(.d)
  .d_diff <- survival::survdiff(survival::Surv(PFI.time.1, PFI.1) ~ group, data = .data)
  # print(.d_diff)
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
  # print(kmp)
  coxp <-  broom::tidy(survival::coxph(survival::Surv(PFI.time.1, PFI.1) ~ exp, data = .data, na.action = na.exclude))
  # print(coxp)
  
  .fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = .data, na.action = na.exclude) 
  # print(.fit)
  # fn_sur_draw(.gene,kmp,.fit)
  data.frame(KMP=kmp,Coxp=coxp) -> .out
  return(.out)
}

### miRNA OS ----
miRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-name) -> miRNA_clinical
miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(survival=purrr::map2(name,data,.up=75,.low=25,fn_survival)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> miRNA_survival_p_OS
readr::write_tsv(miRNA_survival_p_OS,file.path(survival_path,"resulttable","OS_miRNA_survival_pvalue.txt"))

for(i in 1:nrow(miRNA_survival_p_OS)){
  .up=75
  .low=25
  gene=miRNA_survival_p_OS$name[i]
  miRNA_clinical %>%
    dplyr::filter(name==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=miRNA_survival_p_OS$KMP[i]
  coxp=miRNA_survival_p_OS$Coxp.p.value[i]
  fig_name <- paste("OS",gene,KMP,"pdf", sep = ".")
  # fig_name <- paste("PFI",gene,KMP,"png", sep = ".")
  fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = draw_for_sur, na.action = na.exclude)
  # fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = .data, na.action = na.exclude) 
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2),"COX p =",signif(coxp, 2)),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend= c(0.8,0.8),
                        legend.labs = c(paste("Low", .low,"%",", n = ",low_n,sep=""),
                                        paste("Up", .up,"%",", n = ",high_n,sep="")),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw()
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path,"Figure5C.survival/OS"), width = 4, height = 4)
}
### miRNA PFI ----
miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(survival=purrr::map2(name,data,.up=75,.low=25,fn_survival_PFI)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> miRNA_survival_p_PFI
readr::write_tsv(miRNA_survival_p_PFI,file.path(survival_path,"resulttable","PFI_miRNA_survival_pvalue.txt"))

for(i in 1:nrow(miRNA_survival_p_PFI)){
  .up=75
  .low=25
  gene=miRNA_survival_p_PFI$name[i]
  miRNA_clinical %>%
    dplyr::filter(name==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=miRNA_survival_p_PFI$KMP[i]
  coxp=miRNA_survival_p_PFI$Coxp.p.value[i]
  # fig_name <- paste("OS",gene,KMP,"png", sep = ".")
  fig_name <- paste("PFI",gene,KMP,"pdf", sep = ".")
  fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = draw_for_sur, na.action = na.exclude)
  # fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = .data, na.action = na.exclude) 
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2),"COX p =",signif(coxp, 2)),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend= c(0.8,0.8),
                        legend.labs = c(paste("Low", .low,"%",", n = ",low_n,sep=""),
                                        paste("Up", .up,"%",", n = ",high_n,sep="")),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw()
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path,"Figure5C.survival/PFI"), width = 4, height = 4)
}
#### mRNA ----
mRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-symbol) -> mRNA_clinical

### for OS ----
mRNA_clinical %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(survival=purrr::map2(symbol,data,.up=75,.low=25,fn_survival)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> mRNA_survival_p_OS

readr::write_tsv(mRNA_survival_p_OS,file.path(survival_path,"resulttable","OS_mRNA_survival_pvalue.txt"))
### OS pic -----
for(i in 1:nrow(mRNA_survival_p_OS)){
  .up=75
  .low=25
  gene=mRNA_survival_p_OS$symbol[i]
  mRNA_clinical %>%
    dplyr::filter(symbol==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=mRNA_survival_p_OS$KMP[i]
  coxp=mRNA_survival_p_OS$Coxp.p.value[i]
  fig_name <- paste("OS",gene,KMP,"pdf", sep = ".")
  # fig_name <- paste("PFI",gene,KMP,"png", sep = ".")
  
  fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = draw_for_sur, na.action = na.exclude)
  # fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = draw_for_sur, na.action = na.exclude)
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2),"Cox p =",signif(coxp, 2)),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend.labs = c(paste("Low", .low,"%",", n = ",low_n,sep=""),
                                        paste("Up", .up,"%",", n = ",high_n,sep="")),
                        legend= c(0.8,0.8),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw()
  ) 
  # dev.off()
  # ggsave(filename = fig_name, device = "png", path = file.path(survival_path,"Figure5C.survival/PFI"), width = 6, height = 6)
  ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path,"Figure5C.survival/OS"), width = 4, height = 4)
}

### PFI -----
mRNA_clinical %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(survival=purrr::map2(symbol,data,.up=75,.low=25,fn_survival_PFI)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> mRNA_survival_p_PFI
readr::write_tsv(mRNA_survival_p_PFI,file.path(survival_path,"resulttable","PFI_mRNA_survival_pvalue.txt"))

### PFI pic -----
for(i in 1:nrow(mRNA_survival_p_PFI)){
  .up=75
  .low=25
  gene=mRNA_survival_p_PFI$symbol[i]
  mRNA_clinical %>%
    dplyr::filter(symbol==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=mRNA_survival_p_PFI$KMP[i]
  coxp=mRNA_survival_p_PFI$Coxp.p.value[i]
  # fig_name <- paste(gene,KMP,"png", sep = ".")
  fig_name <- paste("PFI",gene,KMP,"pdf", sep = ".")
  
  # fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = draw_for_sur, na.action = na.exclude)
  fit <- survival::survfit(survival::Surv(PFI.time.1, PFI.1) ~ group, data = draw_for_sur, na.action = na.exclude)
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2),"Cox p =",signif(coxp, 2)),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend.labs = c(paste("Low", .low,"%",", n = ",low_n,sep=""),
                                        paste("Up", .up,"%",", n = ",high_n,sep="")),
                        legend= c(0.8,0.8),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw()
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path(survival_path,"Figure5C.survival/PFI"), width = 6, height = 6)
  # ggsave(filename = fig_name, device = "png", path = file.path(survival_path,"Figure5C.survival/OS"), width = 6, height = 6)
}


# draw pic -------
mRNA_survival_p %>%
  dplyr::filter(KMP<=0.05) %>%
  dplyr::mutate(logP=-log10(KMP)) %>%
  dplyr::arrange(Worse,logP) %>%
  ggplot(aes(x=Coxp.term,y=symbol,color=Worse,size=logP)) +
  geom_point() +
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
  )  ->p;p
ggsave(file.path(survival_path,"all_gene_survival-point.pdf"),p,width = 8,height = 4)
miRNA_survival_p %>%
  dplyr::filter(KMP<=0.05) %>%
  dplyr::mutate(logP=-log10(KMP)) %>%
  dplyr::arrange(Worse,logP) %>%
  dplyr::rename("miRNA"="name") %>%
  ggplot(aes(x=Coxp.term,y=miRNA,color=Worse,size=logP)) +
  geom_point() +
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
  )  ->p;p
ggsave(file.path(survival_path,"all_miRNA_survival-point.pdf"),p,width = 6,height = 4)


# 联合survival --------------------------------------------------------------
# #### gene-gene -------------------
data <- mRNA_clinical 
  .up=50
  .low=50
  .gene1 = "EZH2"
  .gene2 = "CBX2"
  .gene3 = "CLDN11"
  data %>%
    dplyr::filter(symbol %in% .gene1) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group1=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group1=ifelse(exp<quantile(exp,probs =.low/100),"Low",group1)) -> .data1

  data %>%
    dplyr::filter(symbol %in% .gene2) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group2=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group2=ifelse(exp<quantile(exp,probs =.low/100),"Low",group2)) -> .data2
  data %>%
    dplyr::filter(symbol %in% .gene3) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group3=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group3=ifelse(exp<quantile(exp,probs =.low/100),"Low",group3)) -> .data3
  .data2 %>%
    dplyr::inner_join(.data1,by="sample") %>%
    dplyr::inner_join(.data3,by="sample") %>%
    dplyr::mutate(group=ifelse(group1=="Up"&group2=="Up"&group3=="Up","Up","NA")) %>%
    dplyr::mutate(group=ifelse(group1=="Low"&group2=="Low"&group3=="Low","Low",group)) %>%
    dplyr::filter(! group=="NA")-> .data
  
  .data2 %>%
    dplyr::inner_join(.data1,by="sample") %>%
    dplyr::inner_join(.data3,by="sample") %>%
    dplyr::mutate(group=ifelse(group1=="Up"&group2=="Up"&group3=="Low","Up","NA")) %>%
    dplyr::mutate(group=ifelse(group1=="Low"&group2=="Low"&group3=="Up","Low",group)) %>%
    dplyr::filter(! group=="NA")-> .data
  ##### PFI survival ----
  .d_diff <- survival::survdiff(survival::Surv(PFI.time.1.x, PFI.1.x) ~ group, data = .data)
  # print(.d_diff)
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
  # print(kmp)
  # coxp <-  broom::tidy(survival::coxph(survival::Surv(PFI.time.1.x, PFI.1.x) ~ exp, data = .data, na.action = na.exclude))
  # print(coxp)
  low_n=.data %>% dplyr::filter(group==paste("Low")) %>% nrow()
  high_n=.data %>% dplyr::filter(group==paste("Up")) %>% nrow()
  .fit <- survival::survfit(survival::Surv(PFI.time.1.x, PFI.1.x) ~ group, data = .data, na.action = na.exclude) 
  # print(.fit)
  # fn_sur_draw(.gene,kmp,.fit)
  fig_name <- paste(.gene1,"+",.gene2,"+",.gene3,"-","PFI",".pdf",sep="")
  library(ggplot2)
  survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(.gene1,.gene2,.gene3, sep = "-"), "Logrank P =", signif(kmp, 2),", PFI"),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend.labs = c(paste(.gene1,",",.gene2," ","Up","+",.gene3," ","Low", .low,"%",", n = ",low_n,sep=""),
                                        paste(.gene1,",",.gene2," ","Low","+",.gene3," ","Up", .up,"%",", n = ",high_n,sep="")),
                        legend= c(0.8,0.8),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw(),
                        font.main = c(16),
                        font.x = c(14),
                        font.y = c(14),
                        font.tickslab = c(12)
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure3"), width = 6, height = 5)
  
  #### OS survival -----
  .d_diff <- survival::survdiff(survival::Surv(OS.x, status.x) ~ group, data = .data)
  # print(.d_diff)
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
  # print(kmp)
  # coxp <-  broom::tidy(survival::coxph(survival::Surv(PFI.time.1.x, PFI.1.x) ~ exp, data = .data, na.action = na.exclude))
  # print(coxp)
  low_n=.data %>% dplyr::filter(group==paste("Low")) %>% nrow()
  high_n=.data %>% dplyr::filter(group==paste("Up")) %>% nrow()
  .fit <- survival::survfit(survival::Surv(OS.x, status.x) ~ group, data = .data, na.action = na.exclude) 
  # print(.fit)
  # fn_sur_draw(.gene,kmp,.fit)
  fig_name <- paste(.gene1,"-",.gene2,"-","OS",".pdf",sep="")
  library(ggplot2)
  survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(.gene1,.gene2, sep = "-"), "Logrank P =", signif(kmp, 2),", OS"),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend.labs = c(paste(.gene1,"-",.gene2," ","Low", .low,"%",", n = ",low_n,sep=""),
                                        paste(.gene1,"-",.gene2," ","Up", .up,"%",", n = ",high_n,sep="")),
                        legend= c(0.8,0.8),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw(),
                        font.main = c(16),
                        font.x = c(14),
                        font.y = c(14),
                        font.tickslab = c(12)
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2"), width = 6, height = 5)
  
  # #### gene-miRNA -------------------
  data1 <- mRNA_clinical 
  data2 <- miRNA_clinical
  .up=50
  .low=50
  .gene5 = "CBX2"
  .gene1 = "EZH2"
  .gene2 = "hsa-miR-101-3p"
  .gene3 = "hsa-miR-30d-5p"
  .gene4 = "hsa-let-7c-5p"
  data1 %>%
    dplyr::filter(symbol %in% .gene5) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group5=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group5=ifelse(exp<quantile(exp,probs =.low/100),"Low",group5)) -> .data5
  
  data1 %>%
    dplyr::filter(symbol %in% .gene1) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group1=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group1=ifelse(exp<quantile(exp,probs =.low/100),"Low",group1)) -> .data1
  
  data2 %>%
    dplyr::filter(name %in% .gene2) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group2=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group2=ifelse(exp<quantile(exp,probs =.low/100),"Low",group2)) -> .data2
  
  data2 %>%
    dplyr::filter(name %in% .gene3) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group3=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group3=ifelse(exp<quantile(exp,probs =.low/100),"Low",group3)) -> .data3
  
  data2 %>%
    dplyr::filter(name %in% .gene4) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group4=ifelse(exp>quantile(exp,probs =.up/100),"Up",NA)) %>%
    dplyr::mutate(group4=ifelse(exp<quantile(exp,probs =.low/100),"Low",group4)) -> .data4
  
  .data2 %>%
    dplyr::inner_join(.data1,by="sample") %>%
    dplyr::inner_join(.data3,by="sample") %>%
    dplyr::inner_join(.data4,by="sample") %>%
    dplyr::inner_join(.data5,by="sample") %>%
    dplyr::mutate(group=ifelse(group1=="Up"& group5=="Up"&group2=="Low" & group3=="Low" & group4=="Low","Up","NA")) %>%
    dplyr::mutate(group=ifelse(group1=="Low"& group5=="Low"&group2=="Up" & group3=="Up" & group4=="Up","Low",group)) %>%
    dplyr::filter(! group=="NA")-> .data
  
  ##### PFI survival ----
  .d_diff <- survival::survdiff(survival::Surv(PFI.time.1.x, PFI.1.x) ~ group, data = .data)
  # print(.d_diff)
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
  # print(kmp)
  # coxp <-  broom::tidy(survival::coxph(survival::Surv(PFI.time.1.x, PFI.1.x) ~ exp, data = .data, na.action = na.exclude))
  # print(coxp)
  low_n=.data %>% dplyr::filter(group==paste("Low")) %>% nrow()
  high_n=.data %>% dplyr::filter(group==paste("Up")) %>% nrow()
  .fit <- survival::survfit(survival::Surv(PFI.time.1.x, PFI.1.x) ~ group, data = .data, na.action = na.exclude) 
  # print(.fit)
  # fn_sur_draw(.gene,kmp,.fit)
  fig_name <- paste(.gene1,"+",.gene5,"+",.gene2,"+",.gene3,"+",.gene4,"-","PFI",".pdf",sep="")
  library(ggplot2)
  survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(.gene1,"+",.gene5,"+",.gene2,"+",.gene3,"+",.gene4, sep = ""), "Logrank P =", signif(kmp, 2),", PFI"),
                        xlab = "Survival in days",
                        ylab = 'Probability of survival',
                        legend.title = "Expression group:",
                        legend.labs = c(paste(.gene1,",",.gene5,"-","Low"," ", .low,"%"," + ",.gene2,",",.gene3,",",.gene4,"-","Up"," ", .low,"%",", n = ",low_n,sep=""),
                                        paste(.gene1,",",.gene5,"-","Up"," ", .up,"%"," + ",.gene2,",",.gene3,",",.gene4,"-","Up"," ", .up,"%",", n = ",high_n,sep="")),
                        legend= c(0.8,0.8),
                        # risk.table = TRUE,
                        # tables.height = 0.2,
                        palette = c("#E7B800", "#2E9FDF"),
                        ggtheme = theme_bw(),
                        font.main = c(16),
                        font.x = c(14),
                        font.y = c(14),
                        font.tickslab = c(12)
  ) 
  # dev.off()
  ggsave(filename = fig_name, device = "pdf", path = file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure3"), width = 6, height = 5)  
# stage -------------------------------------------------------------------
fn_kruskal <- function(.data){
  .data %>%
    tidyr::drop_na() %>%
    dplyr::mutate(stage=as.factor(stage))-> .data
  broom::tidy(kruskal.test(exp~stage,data=.data)) -> .out
}
mRNA_clinical %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(kruskal=purrr::map(data,fn_kruskal)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> mRNA_stage_p
mRNA_stage_p %>%
  dplyr::filter(p.value<=0.05) -> mRNA_stage_p.sig
readr::write_tsv(mRNA_stage_p,file.path(survival_path,"resulttable","mRNA_stage_pvalue.txt"))

miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(kruskal=purrr::map(data,fn_kruskal)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> miRNA_stage_p
readr::write_tsv(miRNA_stage_p,file.path(survival_path,"resulttable","miRNA_stage_pvalue.txt"))

# draw pic ----
comp_list <- list(c("Normal(TA)","stage I"),c("stage I", "stage II"), c("stage II", "stage III"),c("stage III", "stage IV"))
miRNA_stage_p %>%
  dplyr::filter(p.value<=0.06) -> miRNA_stage_p.sig
miRNA_clinical %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::select(name,sample,exp,stage) %>%
  rbind(miRNA_exp_normal.gather) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  dplyr::filter(name %in% miRNA_stage_p.sig$name) %>%
  dplyr::filter(! name %in% EZH2upstream_miRNA_list) %>%
  dplyr::arrange(stage) %>%
  ggpubr::ggboxplot(x = "stage", y = "log2exp",
                    color = "stage", palette = "npg", add = "jitter",
                    facet.by = "name") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 19) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(15, 16,17,18))->p;p
ggsave(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4/Figure5","miRNA_stage-box.pdf"),p,width = 10,height = 3)


mRNA_exp_normal.gather %>%
  dplyr::filter(symbol %in% mRNA_stage_p.sig$symbol) -> mRNA_exp_normal.gather.filter
mRNA_clinical %>%
  dplyr::filter(symbol %in% mRNA_stage_p.sig$symbol) %>%
  dplyr::filter(! symbol %in% keys) %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::select(symbol,sample,exp,stage) %>%
  rbind(mRNA_exp_normal.gather.filter) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  dplyr::arrange(stage) %>%
  ggpubr::ggboxplot(x = "stage", y = "log2exp",
                    color = "stage", palette = "npg", add = "jitter",
                    facet.by = "symbol") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 20) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(16, 16.5,17,17.5))->p;p
ggsave(file.path(survival_path,"all_mRNA_stage-box.pdf"),p,width = 20,height = 15)


# metastic ----------------------------------------------------------------
fn_wilcoxon <- function(.data){
  .data %>%
    dplyr::select(sample,exp,metastasis) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(stage=as.factor(metastasis))-> .data
  broom::tidy(wilcox.test(exp~metastasis,data=.data)) -> .out
}
mRNA_clinical %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(kruskal=purrr::map(data,fn_wilcoxon)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> mRNA_metastic_p
mRNA_metastic_p %>%
  dplyr::filter(p.value<=0.05) -> mRNA_metastic_p.sig
readr::write_tsv(mRNA_metastic_p,file.path(survival_path,"resulttable","mRNA_metastic_pvalue.txt"))

miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(kruskal=purrr::map(data,fn_wilcoxon)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> miRNA_metastic_p
readr::write_tsv(miRNA_metastic_p,file.path(survival_path,"resulttable","miRNA_metastic_pvalue.txt"))

# draw pic ----
comp_list_metas <- list(c("Normal(TA)","m0"),c("m0", "m1"))
miRNA_clinical %>%
  tidyr::unnest() %>%
  dplyr::select(name,sample,exp,metastasis) %>%
  tidyr::drop_na() -> miRNA_clinical.filter
miRNA_exp_normal.gather %>% dplyr::rename("metastasis"="stage") %>%
  rbind(miRNA_clinical.filter) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  ggpubr::ggboxplot(x = "metastasis", y = "log2exp",
                    color = "metastasis", palette = "npg", add = "jitter",
                    facet.by = "name") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 19) +
  ggpubr::stat_compare_means(comparisons = comp_list_metas,method = "wilcox.test",label.y = c(16,17,18))->p;p
ggsave(file.path(survival_path,"all_miRNA_metastic-box.pdf"),p,width = 10,height = 10)



mRNA_exp_normal.gather %>%
  dplyr::rename("metastasis"="stage") -> mRNA_exp_normal.gather.filter
mRNA_clinical %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::filter(metastasis=="m0") %>%
  dplyr::select(symbol,sample,exp,metastasis) -> mRNA_clinical.metas.m0
mRNA_clinical %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::filter(metastasis=="m1") %>%
  dplyr::select(symbol,sample,exp,metastasis) -> mRNA_clinical.metas.m1

mRNA_exp_normal.gather.filter %>%
  rbind(mRNA_clinical.metas.m0) %>%
  rbind(mRNA_clinical.metas.m1) %>%
  dplyr::filter(symbol %in% mRNA_metastic_p.sig$symbol) %>%
  dplyr::filter(! symbol %in% keys) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  ggpubr::ggboxplot(x = "metastasis", y = "log2exp",
                    color = "metastasis", palette = "npg", add = "jitter",
                    facet.by = "symbol") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 20) +
  ggpubr::stat_compare_means(comparisons = comp_list_metas,method = "wilcox.test",label.y = c(16,17,18))->p;p
ggsave(file.path(survival_path,"all_mRNA_metastic-box.pdf"),p,width = 12,height = 10)

#### for EHZ2 and CBX2 ----
mRNA_exp_normal.gather.filter %>%
  rbind(mRNA_clinical.metas.m0) %>%
  rbind(mRNA_clinical.metas.m1) %>%
  dplyr::filter(symbol %in% c("EZH2","CBX2")) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  ggpubr::ggboxplot(x = "metastasis", y = "log2exp",
                    color = "metastasis", palette = "npg", add = "jitter",
                    facet.by = "symbol") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 20) +
  ggpubr::stat_compare_means(comparisons = comp_list_metas,method = "wilcox.test",label.y = c(16,17,18))->p;p
ggsave(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2-CBX2_metastic-box.pdf"),p,width = 5,height = 3)

#### for E2F1 and SOX4 ----
mRNA_exp_normal.gather.filter %>%
  rbind(mRNA_clinical.metas.m0) %>%
  rbind(mRNA_clinical.metas.m1) %>%
  dplyr::filter(symbol %in% c("E2F1","SOX4")) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  ggpubr::ggboxplot(x = "metastasis", y = "log2exp",
                    color = "metastasis", palette = "npg", add = "jitter",
                    facet.by = "symbol") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 20) +
  ggpubr::stat_compare_means(comparisons = comp_list_metas,method = "wilcox.test",label.y = c(16,17,18))->p;p
ggsave(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2","EZH2-CBX2_metastic-box.pdf"),p,width = 5,height = 3)


# info combination --------------------------------------------------------

miRNA_survival_p_OS %>%
  dplyr::select(name,KMP,Coxp.p.value,Worse) %>%
  dplyr::rename("OS_KMP"="KMP","OS_Cox"="Coxp.p.value") %>%
  tidyr::gather(-name,-Worse,key="Term",value="p.value") -> miRNA_survival_p_os.gather
miRNA_survival_p_PFI %>%
  dplyr::select(name,KMP,Coxp.p.value,Worse) %>%
  dplyr::rename("PFI_KMP"="KMP","PFI_Cox"="Coxp.p.value") %>%
  tidyr::gather(-name,-Worse,key="Term",value="p.value") -> miRNA_survival_p_PFI.gather
miRNA_stage_p %>%
  dplyr::mutate(Term="Stage",Worse="") %>%
  dplyr::select(name,Worse,Term,p.value) %>%
  dplyr::ungroup() -> miRNA_stage_p.gather
miRNA_metastic_p %>%
  dplyr::mutate(Term="Metastasis",Worse="") %>%
  dplyr::select(name,Worse,Term,p.value) %>%
  dplyr::ungroup() -> miRNA_metatic_p.gather
miRNA_survival_p_os.gather %>%
  rbind(miRNA_survival_p_PFI.gather) %>%
  rbind(miRNA_stage_p.gather) %>%
  rbind(miRNA_metatic_p.gather) -> miRNA_clinical_result

miRNA_clinical_result %>%
  dplyr::filter(! name %in% EZH2upstream_miRNA_list) %>%
  dplyr::mutate(`-log10P`=-log10(p.value)) %>%
  dplyr::filter(p.value<=0.05) -> tmp
tmp %>%
  ggplot(aes(x=Term,y=name,color=Worse)) +
  geom_tile(aes(fill = `-log10P`,colour = Worse)) +
  scale_fill_gradientn(colours=c("#00BFFF","red"),
                       limits=c(1,4),
                       name = "-log10(P)") +
  scale_x_discrete(limits=c("Stage","OS_KMP","PFI_KMP","PFI_Cox")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure4/Figure5","miRNAs_clinical.pdf"),width=5,height=4)


mRNA_survival_p_OS %>%
  dplyr::select(symbol,KMP,Coxp.p.value,Worse) %>%
  dplyr::rename("OS_KMP"="KMP","OS_Cox"="Coxp.p.value") %>%
  tidyr::gather(-symbol,-Worse,key="Term",value="p.value") -> mRNA_survival_p_OS.gather
mRNA_survival_p_PFI %>%
  dplyr::select(symbol,KMP,Coxp.p.value,Worse) %>%
  dplyr::rename("PFI_KMP"="KMP","PFI_Cox"="Coxp.p.value") %>%
  tidyr::gather(-symbol,-Worse,key="Term",value="p.value") -> mRNA_survival_p_PFI.gather
mRNA_stage_p %>%
  dplyr::mutate(Term="Stage",Worse="") %>%
  dplyr::select(symbol,p.value,Worse,Term) %>%
  dplyr::ungroup() -> mRNA_stage_p.gather
mRNA_metastic_p %>%
  dplyr::mutate(Term="Metastasis",Worse="") %>%
  dplyr::select(symbol,p.value,Worse,Term) %>%
  dplyr::ungroup() -> mRNA_metastic_p.gather
mRNA_survival_p_OS.gather %>%
  rbind(mRNA_survival_p_PFI.gather) %>%
  rbind(mRNA_stage_p.gather) %>%
  rbind(mRNA_metastic_p.gather) -> mRNA_clinical_result

mRNA_clinical_result %>%
  dplyr::mutate(`-log10P`=-log10(p.value)) %>%
  dplyr::filter(p.value<=0.05) -> tmp
tmp %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(symbol,n) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::arrange(n) -> generank.tmp

tmp %>%
  ggplot(aes(x=Term,y=symbol,color=Worse)) +
  geom_point(aes(size = `-log10P`,colour = Worse)) +
  scale_color_manual(values=c("#3CB371", "#CD0000", "#4F94CD")) +
  scale_size_continuous(limits=c(1,5),
                       name = "-log10(P)") +
  scale_x_discrete(limits=c("PFI_Cox","PFI_KMP","OS_Cox","OS_KMP","Metastasis","Stage")) +
  scale_y_discrete(limits=generank.tmp$symbol) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2)
  )
ggsave(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure5","mRNAs_clinical.pdf"),width=6,height=10)
