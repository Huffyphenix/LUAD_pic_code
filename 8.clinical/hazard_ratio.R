# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

clinical_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/鐢熷瓨鍒嗘瀽/data/LUAD"
clinical_path_1 <- "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/survival"
survival_path <- "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/Figure/Figure5"
data_path<- "H:/data"
exp_path <- "H:/data/TCGA/TCGA_data"
enrich_path <- "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/閫氳矾瀵岄泦/LUAD-noFC-prob0.9-kegg-gsea"
target_path <- "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"


# load data ---------------------------------------------------------------

clinical <- readr::read_tsv(file.path(clinical_path,"clinical_info_multi_survival.txt"))
survival <- readr::read_tsv(file.path(clinical_path_1,"cancer_cell_survival_time_table.txt"))
survival %>%
  dplyr::mutate(sample=substr(bcr_patient_barcode,9,12)) %>%
  dplyr::select(sample,PFI.1,PFI.time.1) %>%
  dplyr::mutate(PFI.time.1=as.numeric(PFI.time.1)/30) -> PFI_survival_time
miRNA_exp <- readr::read_rds(file.path(exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()

targets <- readr::read_tsv(file.path(target_path,"all_EHZ2_CBX2_common_targets.DE_info"))
targets %>%
  dplyr::filter(prob>0.99 & log2FC <= -0.585) %>%
  dplyr::filter(case_mean>=30) -> targets_down
ppar <- readr::read_tsv(file.path(enrich_path,"gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description == "PPAR signaling pathway")%>%
  .$SYMBOL
keys <- c("EZH2","CBX2","E2F1","SOX4","JUN")

# filter ------------------------------------------------------------------

gene_exp %>%
  dplyr::filter(symbol %in% c(ppar)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==0) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() -> mRNA_exp.gather
gene_exp %>%
  dplyr::filter(symbol %in% c(ppar)) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="exp") %>%
  dplyr::filter(substr(sample,14,14)==1) %>%
  dplyr::mutate(sample=substr(sample,9,12)) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(stage="Normal(TA)") -> mRNA_exp_normal.gather

clinical %>%
  dplyr::mutate(status=ifelse(vital_status=="dead",as.integer(1),as.integer(0))) %>%
  dplyr::mutate(OS=ifelse(is.na(days_to_death),days_to_last_followup,days_to_death)) %>%
  dplyr::select(patient_id,status,OS) %>%
  dplyr::rename("sample"="patient_id")-> clinical_data

clinical_data %>% 
  dplyr::full_join(PFI_survival_time,by="sample") -> clinical_data

mRNA_exp.gather %>%
  dplyr::inner_join(clinical_data,by="sample") %>%
  tidyr::nest(-symbol) -> mRNA_clinical


# calculte hazrd ratio--------------------------------------------------------------

mRNA_clinical %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = symbol,
    .f = function(.x, .y){
      print(.y)
      .x %>% 
        dplyr::mutate(group = dplyr::case_when(
          exp <= quantile(exp, 0.4) ~ "1L",
          exp > quantile(exp, 0.6) ~ "2H",
          TRUE ~ "M"
        )) %>%
        dplyr::select(-status,-OS) %>%
        dplyr::filter(PFI.time.1<=60) %>%
        tidyr::drop_na() %>% 
        dplyr::filter(group != "M")-> .d
      .d %>% dplyr::filter(group != "M") -> .dd
      
      survival::coxph(survival::Surv(PFI.time.1,PFI.1 ) ~ group, data = .d) -> .cox
      summary(.cox) -> .z
      survival::survdiff(survival::Surv(PFI.time.1, PFI.1) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      tibble::tibble(
        n = .z$n,
        hr = .z$conf.int[1],
        hr_l = .z$conf.int[3],
        hr_h = .z$conf.int[4],
        coxp = .z$waldtest[3],
        kmp = kmp)
      
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest() -> ppar_hazard_ratio;ppar_hazard_ratio

ppar_hazard_ratio %>% 
  dplyr::mutate(coxp = -log10(coxp)) %>%
  dplyr::mutate(symbol = reorder(symbol,hr,sort))-> plot_ready


# draw pic ----------------------------------------------------------------

plot_ready %>% 
  ggplot(aes(y = hr, x = symbol, ymin=hr_l,ymax=hr_h)) +
  geom_pointrange() +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  coord_flip() +
  ggthemes::theme_gdocs() +
  labs(y = "Hazard Ratio", x = "Symbols") -> p;p
ggsave("F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/Figure/Figure1/PPAR_hazard_ratio.pdf"  )
ggsave(filename = "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/Figure/Figure1/PPAR_hazard_ratio.pdf",device="pdf",width=4,height=3)
ggsave(filename = "F:/鎴戠殑鍧氭灉浜?/ENCODE-TCGA-LUAD/Figure/Figure1/PPAR_hazard_ratio.tiff",device="tiff",width=4,height=3)

knitr::kable(plot_ready %>% 
               dplyr::rename(`hazard ratio` = hr) %>% 
               dplyr::mutate(coxp = 10^-coxp) %>% 
               dplyr::arrange(-`hazard ratio`) %>% 
               dplyr::mutate(hr_l = signif(hr_l, 3), hr_h = signif(hr_h,3)) %>% 
               tidyr::unite(`95%CI`, hr_l, hr_h, sep = "-"))
