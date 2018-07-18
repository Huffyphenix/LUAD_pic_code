library(magrittr)
# data path ---------------------------------------------------------------

miRNA_exp_path <- "H:/data/TCGA/TCGA_data"
data_path <- "H:/data/GSCALite/GTEx"
# laod data ---------------------------------------------------------------

miRNA_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_mirna_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
gene_exp <- readr::read_rds(file.path(miRNA_exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest() 
lung_normal <- readr::read_rds(file.path(data_path,"GTEx_lung_mRNASeq.rds")) %>%
  tidyr::unnest()

# filter data -------------------------------------------------------------
gene_list <- c("CBX2","EZH2")
gene_exp %>%
  dplyr::filter(symbol %in% gene_list) %>%
  dplyr::select(-cancer_types,-entrez_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  dplyr::mutate(sample=substr(sample,9,16)) %>%
  dplyr::as_tibble() %>%
  tidyr::spread(key=symbol,value=gene_exp) %>%
  dplyr::mutate(sample_type=ifelse(substr(sample,6,6)==0,"TCGA Tumor","TCGA Normal"))-> EZH2_CBX2_exp.gather
lung_normal %>%
  dplyr::filter(symbol %in% gene_list) %>%
  dplyr::select(-SMTS,-ensembl_gene_id) %>%
  tidyr::gather(-symbol,key="sample",value="gene_exp") %>%
  tidyr::spread(key=symbol,value=gene_exp) %>%
  dplyr::mutate(sample_type="GTEx Normal")-> EZH2_CBX2_exp.gather.gtex
EZH2_CBX2_exp.gather %>%
  dplyr::filter(sample_type=="TCGA Tumor") -> EZH2_CBX2_exp.gather.T
EZH2_CBX2_exp.gather %>%
  dplyr::filter(sample_type=="TCGA Normal") -> EZH2_CBX2_exp.gather.N

rbind(EZH2_CBX2_exp.gather.T,EZH2_CBX2_exp.gather.N) %>%
  rbind(EZH2_CBX2_exp.gather.gtex) -> EZH2_CBX2_exp.gather.allsamples
# calculation -------------------------------------------------------------
human_read <- function(.x){
  if (.x > 0.1) {
    .x %>% signif(digits = 2) %>% toString()
  } else if (.x < 0.1 && .x > 0.001 ) {
    .x %>% signif(digits = 1) %>% toString()
  } else {
    .x %>% format(digits = 2, scientific = TRUE)
  }
}

broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.T$CBX2,EZH2_CBX2_exp.gather.T$EZH2,method = "spearman")
  ) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=8,y=11,sample_type="TCGA Tumor",n=nrow(EZH2_CBX2_exp.gather.T), estimate = signif(estimate,2)) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value, human_read)) %>% 
  dplyr::mutate(label=purrr::map(
    .x = p.value,
    .y = estimate,
    .z = n,
    .f = function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  ) )-> text.T
  
broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.gtex$CBX2,EZH2_CBX2_exp.gather.gtex$EZH2,method = "spearman")
) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=2.7,y=3,sample_type="GTEx Normal", n=nrow(EZH2_CBX2_exp.gather.gtex), estimate = signif(estimate,2)) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value, human_read)) %>% 
  dplyr::mutate(label=purrr::map(
    .x = p.value,
    .y = estimate,
    .z = n,
    .f = function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  ) ) -> text.GN
broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.N$CBX2,EZH2_CBX2_exp.gather.N$EZH2,method = "spearman")
) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=2.7,y=3,sample_type="TCGA Normal", n=nrow(EZH2_CBX2_exp.gather.N), estimate = signif(estimate,2)) %>%
  dplyr::mutate(p.value = purrr::map_chr(p.value, human_read)) %>% 
  dplyr::mutate(label=purrr::map(
    .x = p.value,
    .y = estimate,
    .z = n,
    .f = function(.x,.y,.z){
      if(grepl(pattern = "e",x=.x)){
        sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
        latex2exp::TeX(glue::glue("r = <<.y>>, p = $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$, n = <<.z>>", .open = "<<", .close = ">>"))
      } else {
        latex2exp::TeX(glue::glue("r = {.y}, p = {.x}, n = {.z}"))
      }
    }
  ) ) -> text.TN
text.T %>%
  rbind(text.TN) %>%
  rbind(text.GN) %>%
  dplyr::as.tbl() %>%
  dplyr::arrange(sample_type)->text

broom::tidy(
  lm(log2(EZH2_CBX2_exp.gather.T$EZH2) ~ log2(EZH2_CBX2_exp.gather.T$CBX2)) -> lm.x
)
plot(log2(EZH2_CBX2_exp.gather.T$EZH2),log2(EZH2_CBX2_exp.gather.T$CBX2))
abline(lm.x)

# EZH2 have no correlation with CBX2 in normal tissues (GTEx data)
library(ggplot2)
EZH2_CBX2_exp.gather.allsamples %>%
  ggplot(aes(x=log2(EZH2),y=log2(CBX2))) +
  geom_point(aes(color = sample_type)) +
  geom_smooth(span = 0.8, se = FALSE, fullrange=TRUE, color = "#039BE5") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  labs(
    x = "EZH2 mRNA (log2)",
    y = "CBX2 mRNA (log2)"
  ) +
  scale_color_manual(
    values = c("#00C5CD", "#00868B", "#EE6363"),
    labels = text$label
  ) +
  facet_wrap(~sample_type,scales = "free") -> p;p
ggsave(filename = "EZH2_CBX2_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 8,height = 3)
ggsave(filename = "EZH2_CBX2_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 8,height = 3)
