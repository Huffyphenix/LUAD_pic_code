.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")
library(magrittr,ggplot2)
# data path ---------------------------------------------------------------

# target_path <- "S:??????/?ҵļ?????/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
target_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
exp_path <- "H:/data/TCGA/TCGA_data"
clinical_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/生存分析/data/LUAD"
clinical_path_1 <- "F:/我的坚果云/ENCODE-TCGA-LUAD/survival"
survival_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure5"
# survival_path <- "S:??????/?ҵļ?????/ENCODE-TCGA-LUAD/Figure/Figure5"
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
survival <- readr::read_tsv(file.path(clinical_path_1,"cancer_cell_survival_time_table.txt"))

EZH2upstream_miRNA_list <- c("hsa-miR-101-3p","hsa-miR-30d-5p","hsa-let-7c-5p")
targetsupstream_miRNA_list <- c("hsa-miR-210-3p","hsa-miR-183-5p","hsa-miR-151a-5p","hsa-miR-1307-3p","hsa-miR-93-5p",
                                "hsa-miR-2355-5p","hsa-miR-141-5p","hsa-miR-141-92b-3p","hsa-miR-141-3p")
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
  # fn_sur_draw(.gene,kmp,.fit)
  data.frame(KMP=kmp,Coxp=coxp) -> .out
  return(.out)
}

miRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-name) -> miRNA_clinical
miRNA_clinical %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(survival=purrr::map2(name,data,.up=75,.low=25,fn_survival)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> miRNA_survival_p
readr::write_tsv(miRNA_survival_p,file.path(survival_path,"resulttable","miRNA_survival_pvalue.txt"))
for(i in 1:nrow(miRNA_survival_p)){
  .up=75
  .low=25
  gene=miRNA_survival_p$name[i]
  miRNA_clinical %>%
    dplyr::filter(name==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=miRNA_survival_p$KMP[i]
  fig_name <- paste(gene,KMP,"png", sep = ".")
  fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = draw_for_sur, na.action = na.exclude)
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2)),
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
  ggsave(filename = fig_name, device = "png", path = file.path(survival_path,"survival"), width = 6, height = 6)
}

#### mRNA ----
mRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-symbol) -> mRNA_clinical
mRNA_clinical %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(survival=purrr::map2(symbol,data,.up=75,.low=25,fn_survival)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Worse=ifelse(Coxp.estimate>0,"H", "L")) %>%
  dplyr::ungroup() -> mRNA_survival_p
readr::write_tsv(mRNA_survival_p,file.path(survival_path,"resulttable","mRNA_survival_pvalue.txt"))
for(i in 1:nrow(mRNA_survival_p)){
  .up=75
  .low=25
  gene=mRNA_survival_p$symbol[i]
  mRNA_clinical %>%
    dplyr::filter(symbol==gene) %>%
    tidyr::unnest() %>%
    dplyr::mutate(group=ifelse(exp>quantile(exp,probs =.up/100),paste("Up", .up,"%",sep=""),NA)) %>%
    dplyr::mutate(group=ifelse(exp<quantile(exp,probs =.low/100),paste("Low",.low,"%",sep=""),group)) -> draw_for_sur
  low_n=draw_for_sur %>% dplyr::filter(group==paste("Low",.low,"%",sep="")) %>% nrow()
  high_n=draw_for_sur %>% dplyr::filter(group==paste("Up",.up,"%",sep="")) %>% nrow()
  KMP=mRNA_survival_p$KMP[i]
  fig_name <- paste(gene,KMP,"png", sep = ".")
  fit <- survival::survfit(survival::Surv(OS, status) ~ group, data = draw_for_sur, na.action = na.exclude)
  # pdf(file.path(survival_path,"survival",fig_name),width = 6, height = 6)
  survminer::ggsurvplot(fit,pval=F, #pval.method = T,
                        surv.median.line = "hv",
                        title = paste(paste(gene, sep = "-"), "Logrank P =", signif(KMP, 2)),
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
  ggsave(filename = fig_name, device = "png", path = file.path(survival_path,"survival"), width = 6, height = 6)
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
miRNA_clinical %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::select(name,sample,exp,stage) %>%
  rbind(miRNA_exp_normal.gather) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  dplyr::arrange(stage) %>%
  ggpubr::ggboxplot(x = "stage", y = "log2exp",
                    color = "stage", palette = "npg", add = "jitter",
                    facet.by = "name") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 19) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(15, 16,17,18))->p;p
ggsave(file.path(survival_path,"all_miRNA_stage-box.pdf"),p,width = 15,height = 10)

mRNA_exp_normal.gather %>%
  dplyr::filter(symbol %in% mRNA_stage_p.sig$symbol) -> mRNA_exp_normal.gather.filter
mRNA_clinical %>%
  dplyr::filter(symbol %in% mRNA_stage_p.sig$symbol) %>%
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
  ggpubr::stat_compare_means(label.y = 18) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label.y = c(11, 12,13,14))->p;p
ggsave(file.path(survival_path,"all_mRNA_stage-box.pdf"),p,width = 18,height = 10)


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
  dplyr::filter(symbol %in% mRNA_metastic_p.sig$symbol) %>%
  dplyr::rename("metastasis"="stage") -> mRNA_exp_normal.gather.filter
mRNA_clinical %>%
  dplyr::filter(symbol %in% mRNA_metastic_p.sig$symbol) %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::filter(metastasis=="m0") %>%
  dplyr::select(symbol,sample,exp,metastasis) -> mRNA_clinical.metas.m0
mRNA_clinical %>%
  dplyr::filter(symbol %in% mRNA_metastic_p.sig$symbol) %>%
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::filter(metastasis=="m1") %>%
  dplyr::select(symbol,sample,exp,metastasis) -> mRNA_clinical.metas.m1

mRNA_exp_normal.gather.filter %>%
  rbind(mRNA_clinical.metas.m0) %>%
  rbind(mRNA_clinical.metas.m1) %>%
  dplyr::mutate(log2exp=log2(exp)) %>%
  ggpubr::ggboxplot(x = "metastasis", y = "log2exp",
                    color = "metastasis", palette = "npg", add = "jitter",
                    facet.by = "symbol") +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(label.y = 20) +
  ggpubr::stat_compare_means(comparisons = comp_list_metas,method = "wilcox.test",label.y = c(17,18,19))->p;p
ggsave(file.path(survival_path,"all_mRNA_metastic-box.pdf"),p,width = 10,height = 10)

