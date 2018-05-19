data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"
genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
genelist <- readr::read_tsv(file.path(genelist_path,"all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(prob>=0.99 & abs(log2FC)>=0.585) %>%
  dplyr::filter(log2FC<0)

common_gene_cnv <- readr::read_tsv(file.path(data_path,"Figure4/Figure5","cna_of_common_targets.txt"),col_names = F) %>%
  t()
colnames(common_gene_cnv) <- common_gene_cnv["X2",]
common_gene_cnv <- common_gene_cnv[-c(1:2),]

common_gene_cnv %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  tidyr::gather(-SAMPLE_ID,key="sample",value="CNV") %>%
  dplyr::mutate(SAMPLE_ID=as.character(SAMPLE_ID)) %>%
  dplyr::mutate(CNV=as.numeric(CNV)) -> common_gene_cnv
common_gene_cnv %>%
  dplyr::filter(SAMPLE_ID %in% genelist$gene_id.x) %>%
  dplyr::filter(CNV==2) %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::select(-sample) %>%
  dplyr::mutate(homo_amp=sum(CNV)/2) %>%
  dplyr::select(-CNV) %>%
  unique() -> common_gene_cnv.homo_amp
common_gene_cnv %>%
  dplyr::filter(SAMPLE_ID %in% genelist$gene_id.x) %>%
  dplyr::filter(CNV==-2) %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::select(-sample) %>%
  dplyr::mutate(homo_del=sum(CNV)/-2) %>%
  dplyr::select(-CNV) %>%
  unique() -> common_gene_cnv.homo_dele
common_gene_cnv$sample %>% unique() %>% length() -> smaple_n

common_gene_cnv.homo_amp %>%
  dplyr::full_join(common_gene_cnv.homo_dele,by="SAMPLE_ID") %>%
  dplyr::mutate(n=smaple_n) %>%
  dplyr::mutate(homo_amp=ifelse(!is.na(homo_amp),100*homo_amp/n,0)) %>%
  dplyr::mutate(homo_del=ifelse(!is.na(homo_del),100*homo_del/n,0)) %>%
  dplyr::arrange(desc(homo_del)) -> DOWN_commone_targets.CNV_percent
DOWN_commone_targets.CNV_percent %>%
  readr::write_tsv(file.path(data_path,"Figure4/Figure5","Down_common_targets.CNV-percent.tsv"))

# draw figure -------------------------------------------------------------

DOWN_commone_targets.CNV_percent %>%
  dplyr::filter(homo_del>5) %>%
  ggplot(aes(x=SAMPLE_ID,y=homo_del,fill=SAMPLE_ID)) +
  geom_bar(stat = "identity") +
  xlab("Gene") +
  ylab("CNV Deletion (%)") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent', color='black'),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15)
  )
ggsave(file.path(data_path,"Figure4/Figure5","Down_common_targets.CNV-percent(del_5).pdf"))
