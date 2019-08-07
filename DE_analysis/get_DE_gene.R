####################### DE analysis between high and low CBX2, EZH2, both high and low
library(magrittr)
library(ggplot2)
library(ggbeeswarm)

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
# path --------------------------------------------------------------------
# HUST
basic_path <- file.path("S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD")
exp_path <- "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
res_path <- file.path(basic_path,"Figure/DE_between_high_low_CBX2_EZH2")

# Ezhou
basic_path <- file.path("E:/我的坚果云/ENCODE-TCGA-LUAD")
exp_path <- "G:/data/TCGA/TCGA_data"
res_path <- file.path(basic_path,"Figure/DE_between_high_low_CBX2_EZH2")

# load data ---------------------------------------------------------------

gene_exp <- readr::read_rds(file.path(exp_path,"pancan33_expr.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()

# each do one time
key_gene <- "CBX2"
key_gene <- "EZH2"
key_gene <- c("CBX2","EZH2")

# classify samples into two groups ----------------------------------------
if(length(key_gene)>1){
  gene_exp %>%
    dplyr::filter(symbol %in% key_gene) %>%
    tidyr::gather(-cancer_types,-symbol, -entrez_id,key="sample",value="exp") %>%
    dplyr::filter(substr(sample,14,15)=="01") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(group = ifelse(exp>quantile(exp,0.5),"high","low")) %>%
    dplyr::ungroup() %>%
    dplyr::select(symbol,sample,group) %>%
    tidyr::spread(key="symbol",value="group") %>%
    dplyr::filter(CBX2 == EZH2) %>%
    dplyr::mutate(group=CBX2) %>%
    dplyr::select(sample,group) -> group_info
  
  out_path <- file.path(res_path,paste(key_gene,collapse = "_and_"))
}else{
  gene_exp %>%
    dplyr::filter(symbol %in% key_gene) %>%
    tidyr::gather(-cancer_types,-symbol, -entrez_id,key="sample",value="exp") %>%
    dplyr::filter(substr(sample,14,15)=="01") %>%
    dplyr::mutate(group = ifelse(exp>quantile(exp,0.5),"high","low")) %>%
    dplyr::select(sample,group) -> group_info
  
  out_path <- file.path(res_path,paste(key_gene))
}
# compare gene expression itself between high and low group ---------------
gene_exp %>%
  dplyr::filter(symbol %in% key_gene) %>%
  tidyr::gather(-cancer_types,-symbol, -entrez_id,key="sample",value="exp") %>%
  dplyr::inner_join(group_info, by="sample") %>%
  dplyr::mutate("log2 (mRNA Exp.)"= log2(exp+0.01)) %>%
  ggplot(aes(x=symbol,y=`log2 (mRNA Exp.)`)) +
  geom_violin() +
  geom_quasirandom(aes(color = group),size=0.2)  +
  my_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10)
  )
ggsave(file.path(out_path,"high_low_exp_comp.tiff"),device = "tiff",height = 3,width = 4)
# DE analysis -------------------------------------------------------------

fn_DE <- function(gene,group,exp){
  exp %>%
    dplyr::filter(entrez_id %in% gene) %>%
    tidyr::gather(-cancer_types,-symbol, -entrez_id,key="sample",value="exp") %>%
    dplyr::inner_join(group, by="sample") -> .tmp
  
  .tmp %>%
    dplyr::filter(group == "high") %>%
    .$exp %>%
    mean() -> mean_high
  
  .tmp %>%
    dplyr::filter(group == "low") %>%
    .$exp %>%
    mean() -> mean_low
  
  broom::tidy(wilcox.test(exp ~ group, data = .tmp)) %>%
    dplyr::mutate(mean_high = mean_high, mean_low = mean_low, log2FC_HvsL = log2((mean_high+1)/(mean_low+1))) %>%
    dplyr::mutate(DE = ifelse(log2FC_HvsL>0,"Upregulate","Downregulate"))
}

gene <- "100130426"
group <- group_info
exp <- gene_exp
rm(gene,group,exp,.tmp,mean_high,mean_low)

gene_exp %>%
  dplyr::select(entrez_id,symbol) %>%
  dplyr::mutate(res = purrr::map(entrez_id,fn_DE,group=group_info,exp=gene_exp)) %>%
  tidyr::unnest() -> DE_res

DE_res %>%
  readr::write_tsv(file.path(out_path,"DE_gene_res.tsv"))



