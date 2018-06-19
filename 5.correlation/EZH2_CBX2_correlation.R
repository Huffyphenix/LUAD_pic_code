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

# Visualize your data using scatter plots
library("ggpubr")
ggscatter(EZH2_CBX2_exp.gather.T, x = "CBX2", y = "EZH2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")

# Preleminary test to check the test assumptions
shapiro.test(EZH2_CBX2_exp.gather.T$CBX2) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(EZH2_CBX2_exp.gather.T$EZH2) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(EZH2_CBX2_exp.gather.gtex$CBX2) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(EZH2_CBX2_exp.gather.gtex$EZH2) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(EZH2_CBX2_exp.gather.N$EZH2) # p-value < 0.05, don't follow a normal distribution.
shapiro.test(EZH2_CBX2_exp.gather.N$CBX2) # p-value < 0.05, don't follow a normal distribution.

#### So we need to do spearman rank correlation.

# calculation -------------------------------------------------------------

## EZH2 have no correlation with CBX2 in normal tissues (GTEx data)
broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.T$CBX2,EZH2_CBX2_exp.gather.T$EZH2,method = "spearman")
  ) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=8,y=11,sample_type="TCGA Tumor") %>%
  dplyr::mutate(label=paste("Spearman r = ",round(estimate,2)," p = ",signif(fdr,2),"\n","n = ",nrow(EZH2_CBX2_exp.gather.T),
                            "                                    ",sep="")) -> text.T
broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.gtex$CBX2,EZH2_CBX2_exp.gather.gtex$EZH2,method = "spearman")
) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=2.7,y=3,sample_type="GTEx Normal") %>%
  # dplyr::mutate(top = c(-trunc(1/fdr))) %>%
  # dplyr::mutate(bottom=round(fdr*10^(-top),2)) %>%
  # dplyr::mutate(p.label=paste(bottom,"*10",sep="")) %>%
  dplyr::mutate(label=paste("Spearman r = ",round(estimate,4)," p = ",signif(fdr,2),"\n","n = ",nrow(EZH2_CBX2_exp.gather.gtex),
                            "                                  ",sep="")) -> text.N


text.T %>%
  rbind(text.N)  %>%
  rbind(text.N.tcga)->text
broom::tidy(
  lm(log2(EZH2_CBX2_exp.gather.T$EZH2) ~ log2(EZH2_CBX2_exp.gather.T$CBX2)) -> lm.x
)
plot(log2(EZH2_CBX2_exp.gather.T$EZH2),log2(EZH2_CBX2_exp.gather.T$CBX2))
abline(lm.x)

library(ggplot2)
EZH2_CBX2_exp.gather.gtex %>%
  ggplot(aes(x=log2(EZH2),y=log2(CBX2))) +
  geom_point(color = "#00BCD4") +
  geom_smooth(span = 0.8, se = FALSE, fullrange=TRUE, color = "#039BE5",method = "auto") +
  facet_wrap(~sample_type,scales = "free") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
) +
  geom_text(data=text.N,aes(x=x,y=y,label=label),hjust=0.5) +
  labs(
    x = "EZH2 mRNA (log2)",
    y = "CBX2 mRNA (log2)"
  ) -> p;p
ggsave(filename = "EZH2_CBX2_GTEx_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 4,height = 3)
ggsave(filename = "EZH2_CBX2_GTEx_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 4,height = 3)

# EZH2 have correlation with CBX2 in tumor tissues (TCGA data)
EZH2_CBX2_exp.gather.T %>%
  ggplot(aes(x=log2(EZH2),y=log2(CBX2))) +
  geom_point(color = c("#FF6A6A")) +
  geom_smooth(span = 0.8, se = FALSE, fullrange=TRUE, color = "#039BE5",method = "auto") +
  facet_wrap(~sample_type,scales = "free") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  geom_text(data=text.T,aes(x=x,y=y,label=label),hjust=0.5) +
  labs(
    x = "EZH2 mRNA (log2)",
    y = "CBX2 mRNA (log2)"
  ) -> p;p
ggsave(filename = "EZH2_CBX2_mRNA_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 4,height = 3)
ggsave(filename = "EZH2_CBX2_mRNA_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 4,height = 3)

# EZH2 dont have correlation with CBX2 in normal tissues (TCGA data)
broom::tidy(
  cor.test(EZH2_CBX2_exp.gather.N$CBX2,EZH2_CBX2_exp.gather.N$EZH2,method = "spearman")
) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(fdr=p.adjust(p.value,method = "fdr")) %>%
  dplyr::mutate(x=6,y=7.5,sample_type="TCGA Normal") %>%
  # dplyr::mutate(top = c(-trunc(1/fdr))) %>%
  # dplyr::mutate(bottom=round(fdr*10^(-top),2)) %>%
  # dplyr::mutate(p.label=paste(bottom,"*10",sep="")) %>%
  dplyr::mutate(label=paste("Spearman r = ",round(estimate,4)," p = ",signif(fdr,2),"\n","n = ",nrow(EZH2_CBX2_exp.gather.N),
                            "                                    ",sep="")) -> text.N.tcga

EZH2_CBX2_exp.gather.N %>%
  ggplot(aes(x=log2(EZH2),y=log2(CBX2))) +
  geom_point(color = c("#00BCD4")) +
  geom_smooth(span = 0.8, se = FALSE, fullrange=TRUE, color = "#039BE5",method = "auto") +
  facet_wrap(~sample_type,scales = "free") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  geom_text(data=text.N.tcga,aes(x=x,y=y,label=label),hjust=0.5) +
  labs(
    x = "EZH2 mRNA (log2)",
    y = "CBX2 mRNA (log2)"
  ) -> p;p
ggsave(filename = "EZH2_CBX2_tcga_normal_correlation.tiff",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "tiff",width = 4,height = 3)
ggsave(filename = "EZH2_CBX2_tcga_normal_correlation.pdf",path = "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2",device = "pdf",width = 4,height = 3)
