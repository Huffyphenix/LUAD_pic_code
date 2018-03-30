library(magrittr)
# data path ---------------------------------------------------------------

miRNA_list_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data/FC2"
miRNA_exp_path <- "H:/data/TCGA/TCGA_data"
gene_list_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD"
net_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/FC2_De_in_LUAD/FFL/up_E2F1-SOX4-mirna-targets.FFL"

# output data -------------------------------------------------------------
out_path_sup <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
out_path_fig <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"


# load data --------------------------------------------------------------

mirna_target_cor <- readr::read_tsv(file.path(out_path_sup,"Supplementary_Table2.up_FC2_mirna-EZH2_CBX2_targets-spearman.xls"))
targets_methy <- readr::read_tsv(file.path(out_path_sup,"CBX2_EZH2_targets_methy_diff_cor.tsv"))
net_attribute <- readr::read_tsv(file.path(net_path,"attribute-net-all.txt"))
network <- readr::read_tsv(file.path(net_path,"network-uniq.txt"))


# filter ------------------------------------------------------------------
targets_methy %>% 
  dplyr::filter(diff>0 & spm <= -0.4) %>%
  .$symbol -> methy_cor_diff_sig

mirna_target_cor %>%
  dplyr::filter(Cor <= -0.4 & p.value <= 0.05) %>%
  tidyr::unite("mirna_gene",c("mirna_id","gene_id")) -> mirna_gene_unite

network %>%
  dplyr::filter(substr(Source,1,3) =="hsa") %>%
  tidyr::unite("mirna_gene",c("Source","Target")) %>%
  dplyr::inner_join(mirna_gene_unite,by="mirna_gene") %>%
  tidyr::separate("mirna_gene",c("Source","Target"),"_") %>%
  dplyr::select(Source,Target,Regulate_type,Cor) %>%
  dplyr::mutate(Evidence="")-> mirna_gene_cor_sig

network %>%
  dplyr::filter(substr(Source,1,3) !="hsa") %>%
  dplyr::mutate(Cor=0) %>%
  rbind(mirna_gene_cor_sig) -> network.mirna_sig
  
c(network.mirna_sig$Source,network.mirna_sig$Target) %>%
  unique() -> genelist_net

data.frame(Gene=c("Methylation",methy_cor_diff_sig) %>% as.character(),gene_type=c(4,3,3,3,3,3,3) %>% as.integer(),DE_trend=2 %>% as.integer()) %>% 
  as_tibble()-> add_methy_and_methygene_for_attr
data.frame(Source=c("Methylation"),Target=methy_cor_diff_sig,Regulate_type=5 %>% as.integer(),Cor=0 %>% as.integer(),Evidence="") %>%
  as_tibble() -> add_methy_and_methygene_for_net

net_attribute %>%
  dplyr::filter(Gene %in% genelist_net) %>%
  rbind(add_methy_and_methygene_for_attr) %>%
  unique() -> attribute_gene

network.mirna_sig %>%
  rbind(add_methy_and_methygene_for_net) -> network_gene

attribute_gene %>%
  readr::write_tsv(file.path(net_path,"attribute_with_methy-mirna.txt"))
network_gene %>%
  readr::write_tsv(file.path(net_path,"network_with_methy-mirna.txt"))

