library(magrittr)
library(ggplot2)
library(ggpubr)
library(survival)

my_theme <- theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)
# 1. GSVA score correlation with TIL ------
# muli-variable survival analysis -----------------------------------------

basic_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD"
clinical_path <- file.path(basic_path,"survival")
exp_path <- "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
res_path <- file.path(basic_path,"Figure/Figure2/survival")

# load data ---------------------------------------------------------------
# survival and clinical info
clinical_info <- readr::read_tsv(file.path(clinical_path,"LUAD_clinical_info_multi_survival.txt")) %>%
  dplyr::mutate(metastasis = ifelse(pathologic_m %in% c("m1","m1a","m1b"),"yes","1no")) %>%
  dplyr::inner_join(
    tibble::tibble(pathologic_stage=c("stage i","stage ia","stage ib",
                                      "stage ii","stage iia","stage iib",
                                      "stage iiia","stage iiib",
                                      "stage iv" ),
                   stage = c("i","i","i",
                             "ii","ii","ii",
                             "iii","iii",
                             "iv"))
  ) %>%
  dplyr::mutate(smoke_year = ifelse(is.na(number_pack_years_smoked),stopped_smoking_year-year_of_tobacco_smoking_onset,number_pack_years_smoked)) %>%
  dplyr::mutate(smoke_year = ifelse(is.na(stopped_smoking_year) & !is.na(year_of_tobacco_smoking_onset) & is.na(smoke_year),
                                    50,
                                    smoke_year)) %>%
  dplyr::mutate(smoke_year = ifelse(!is.na(stopped_smoking_year) & is.na(year_of_tobacco_smoking_onset) & is.na(smoke_year),
                                    NA,
                                    smoke_year)) %>%
  dplyr::mutate(smoke_year = ifelse(is.na(stopped_smoking_year) & is.na(year_of_tobacco_smoking_onset) & is.na(smoke_year),
                                    0,
                                    smoke_year)) %>%
  dplyr::mutate(smoke_year_group = ifelse(smoke_year>20,">20y",".<=20y")) %>%
  dplyr::mutate(smoke_year_group = ifelse(is.na(smoke_year),NA,smoke_year_group)) %>%
  dplyr::mutate(kras_mut = ifelse(!is.na(kras_mutation_result),"yes","1wild_type")) %>%
  dplyr::mutate(egfr_mut = ifelse(!is.na(egfr_mutation_result),"yes","1wild_type")) %>%
  dplyr::mutate(OS = ifelse(!is.na(days_to_last_followup),days_to_last_followup,days_to_death)) %>%
  dplyr::mutate(OS_status = ifelse(vital_status=="alive",0,1)) %>%
  dplyr::mutate(patient_id = toupper(patient_id)) %>%
  dplyr::select(patient_id,OS,OS_status,metastasis,stage,smoke_year_group,kras_mut,egfr_mut)

PFS_clinical <- readr::read_tsv(file.path(clinical_path,"cancer_cell_survival_time_table.txt"))

PFS_clinical %>%
  dplyr::mutate(patient_id = substr(bcr_patient_barcode,9,12)) %>%
  dplyr::left_join(clinical_info,by="patient_id") %>%
  dplyr::mutate(PFS=as.numeric(PFI.time.1),PFS_status=PFI.1) %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS_status,OS,OS_status,metastasis,stage,smoke_year_group,kras_mut,egfr_mut) -> PFS_OS_clinical

# expression data
gene_exp <- readr::read_rds(file.path(exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() %>%
    dplyr::filter(symbol %in% c("CBX2","EZH2"))

gene_exp %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="barcode",value="exp") %>%
  tidyr::spread(key="symbol",value="exp") %>%
  dplyr::mutate(n=c(1:576)) %>%
  tidyr::gather(-barcode,-n,key="symbol",value="exp") %>%
  dplyr::filter(substr(barcode,14,14)=="0") %>%
  dplyr::mutate(bcr_patient_barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(PFS_OS_clinical,by="bcr_patient_barcode") -> PFS_OS_clinical_exp


# survival analysis -------------------------------------------------------

# 3.2.function to do cox survival analysis  ----
# 3.2.1.univariable cox analysis ---------
fn_survival_test <- function(data,feature){
  print(feature)
  .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
  summary(.cox) -> .z
  
  # KM pvalue
  kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)$chisq,df = length(levels(as.factor(data$group))) - 1)
  
  # mearge results
  tibble::tibble(
    group = rownames(.z$conf.int),
    n = .z$n,
    coef = .z$coefficients[,1],
    hr = .z$conf.int[1],
    hr_l = .z$conf.int[3],
    hr_h = .z$conf.int[4],
    coxp = .z$waldtest[3],
    kmp = kmp) %>%
    dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  # multi-variable cox analysis
}

# 3.2.2.multi-variable cox analysis ---------
fn_survival_test.multiCox <- function(data,uni_sig_feature){
  covariates <- uni_sig_feature$Features %>% unique()
  if(length(covariates)>=1){
    multi_formulas <- as.formula(paste('Surv(time, status)~', paste(covariates,collapse = " + ")))
    model <- coxph(multi_formulas, data = data)
    model.s <- summary(model)
    
    # Extract data
    tibble::tibble(
      Features = rownames(model.s$coefficients),
      n = model.s$n,
      coef = model.s$coefficients[,1],
      hr = model.s$conf.int[,1],
      hr_l = model.s$conf.int[,3],
      hr_h = model.s$conf.int[,4],
      coxp = model.s$coefficients[,5]) %>%
      dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  }else{
    tibble::tibble()
  }
}

# 3.2.3.drawing cox plot ------
fn_cox_plot_1 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    # facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",subtitle = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}

fn_cox_plot_2 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black")
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",title = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}

# CBX2 and EZH2 alone ------
# uni-variable survival analysis

PFS_OS_clinical_exp %>%
  dplyr::filter(PFS <= 30*60) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(group = ifelse(exp > quantile(exp,0.5),"high","1low")) %>%
  dplyr::ungroup() %>%
  dplyr::select(-barcode,-exp) %>%
  tidyr::spread(key="symbol",value="group") %>%
  tidyr::gather(-bcr_patient_barcode,-PFS,-PFS_status,-OS,-OS_status,-n,key="Features",value="group") %>%
  dplyr::mutate(time =PFS, status = PFS_status) %>%
  tidyr::nest(-Features) %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> res.univarite.surv.PFS

res.univarite.surv.PFS %>%
  readr::write_tsv(file.path(res_path,"res.univarite.surv.PFS.tsv"))

# multi-variable survival analysis

res.univarite.surv.PFS %>%
  dplyr::select(Features) %>%
  unique() %>%
  dplyr::mutate(x="all") %>%
  tidyr::nest(-x,.key="sig_features") -> Features_to_do_multi

PFS_OS_clinical_exp %>%
  dplyr::filter(PFS <= 30*60) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(group = ifelse(exp > quantile(exp,0.5),"high","1low")) %>%
  dplyr::ungroup() %>%
  dplyr::select(-barcode,-exp) %>%
  tidyr::spread(key="symbol",value="group") %>%
  dplyr::mutate(time =PFS, status = PFS_status) %>%
  dplyr::mutate(x="all") %>%
  tidyr::nest(-x) %>%
  dplyr::inner_join(Features_to_do_multi,by="x") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(data,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-data,-sig_features) %>%
  tidyr::unnest() -> res.univarite.surv.PFS.multi

res.univarite.surv.PFS.multi %>%
  readr::write_tsv(file.path(res_path,"res.multi-varite.surv.PFS.tsv"))

res.univarite.surv.PFS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-x) %>%
  dplyr::mutate(filename = paste("2.PFS.multi-variable.cox","CBX2_EZH2_sep",sep = "-")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "",title="Multi-variable cox model\nCBX2 and EZH2 (5y PFS)",dir=file.path(res_path),w=6,h=3))

# CBX2 and EZH2 togather ------
# uni-variable survival analysis

PFS_OS_clinical_exp %>%
  dplyr::filter(PFS <= 30*60) %>%
  tidyr::spread(key="symbol",value=exp) %>%
  dplyr::mutate(CBX2 = ifelse(CBX2 > quantile(CBX2,0.5),"high","1low")) %>%
  dplyr::mutate(EZH2 = ifelse(EZH2 > quantile(EZH2,0.5),"high","1low")) %>%
  dplyr::filter(CBX2==EZH2) %>%
  dplyr::mutate(CBX2_EZH2 = paste("CBX2",CBX2,"EZH2",EZH2,sep="_")) %>%
  dplyr::select(-barcode,-CBX2,-EZH2) %>%
  tidyr::gather(-bcr_patient_barcode,-PFS,-PFS_status,-OS,-OS_status,-n,key="Features",value="group") %>%
  dplyr::mutate(time =PFS, status = PFS_status) %>%
  tidyr::nest(-Features) %>%
  dplyr::mutate(surv_res = purrr::map2(data,Features,fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> combine_res.univarite.surv.PFS

combine_res.univarite.surv.PFS %>%
  readr::write_tsv(file.path(res_path,"combine_res.univarite.surv.PFS.tsv"))

# multi-variable survival analysis

combine_res.univarite.surv.PFS %>%
  dplyr::select(Features) %>%
  unique() %>%
  dplyr::mutate(x="all") %>%
  tidyr::nest(-x,.key="sig_features") -> combine_Features_to_do_multi

PFS_OS_clinical_exp %>%
  dplyr::filter(PFS <= 30*60) %>%
  tidyr::spread(key="symbol",value=exp) %>%
  dplyr::mutate(CBX2 = ifelse(CBX2 > quantile(CBX2,0.5),"high","1low")) %>%
  dplyr::mutate(EZH2 = ifelse(EZH2 > quantile(EZH2,0.5),"high","1low")) %>%
  dplyr::filter(CBX2==EZH2) %>%
  dplyr::mutate(CBX2_EZH2 = paste(CBX2,"both",sep = "_")) %>%
  dplyr::select(-barcode,-CBX2,-EZH2) %>%
  dplyr::mutate(time =PFS, status = PFS_status) %>%
  dplyr::mutate(x="all") %>%
  tidyr::nest(-x) %>%
  dplyr::inner_join(combine_Features_to_do_multi,by="x") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(data,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-data,-sig_features) %>%
  tidyr::unnest() -> combine_res.univarite.surv.PFS.multi

combine_res.univarite.surv.PFS.multi %>%
  readr::write_tsv(file.path(res_path,"combine_res.multi-varite.surv.PFS.tsv"))

combine_res.univarite.surv.PFS.multi %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<6) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-6),-6,hr_l),hr_h=ifelse(hr_h>6,6,hr_h)) %>%
  tidyr::nest(-x) %>%
  dplyr::mutate(filename = paste("2.PFS.multi-variable.cox","CBX2_EZH2_combine",sep = "-")) %>%
  dplyr::mutate(res = purrr::map2(data,filename,fn_cox_plot_1,facet = "",title="Multi-variable cox model\nCBX2 and EZH2 (5y PFS)",dir=file.path(res_path),w=6,h=3))
