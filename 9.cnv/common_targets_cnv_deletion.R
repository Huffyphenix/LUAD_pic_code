library(magrittr)
library(ggplot2)
# data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"
data_path<- "H:/data/GSCALite/TCGA/cnv"
result_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/"

# gene list 1 -------------------------------------------------------------

genelist_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new"
genelist <- readr::read_tsv(file.path(genelist_path,"all_EHZ2_CBX2_common_targets.DE_info")) %>%
  dplyr::filter(Class == "Down") %>%
  .$gene_id.x

# gene list 2 -------------------------------------------------------------
enrich_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
cell_cycle_relate <- c("Cell cycle","Oocyte meiosis","DNA replication",
                                   "Homologous recombination","p53 signaling pathway",
                                   "Progesterone-mediated oocyte maturation")
genelist <- readr::read_tsv(file.path(enrich_path,"gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description %in% cell_cycle_relate) %>%
  .$SYMBOL


# gene list 3 -------------------------------------------------------------
enrich_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/通路富集/LUAD-noFC-prob0.9-kegg-gsea"
ppar_relate <- c("PPAR signaling pathway")
genelist <- readr::read_tsv(file.path(enrich_path,"gseaKEGG_result-gather.tsv")) %>%
  dplyr::filter(Description %in% ppar_relate) %>%
  .$SYMBOL


# load cnv data -----------------------------------------------------------

luad_cnv <- readr::read_rds(file.path(data_path,"pancan34_cnv_percent.rds.gz")) %>%
  dplyr::filter(cancer_types=="LUAD") %>%
  tidyr::unnest()

# common_gene_cnv <- readr::read_tsv(file.path(data_path,"Figure4/Figure5","cna_of_common_targets.txt"),col_names = F) %>%
#   t()
# colnames(common_gene_cnv) <- common_gene_cnv["X2",]
# common_gene_cnv <- common_gene_cnv[-c(1:2),]
# 
# common_gene_cnv %>%
#   as.data.frame() %>%
#   dplyr::as.tbl() %>%
#   tidyr::gather(-SAMPLE_ID,key="sample",value="CNV") %>%
#   dplyr::mutate(SAMPLE_ID=as.character(SAMPLE_ID)) %>%
#   dplyr::mutate(CNV=as.numeric(CNV)) -> common_gene_cnv
# common_gene_cnv %>%
#   dplyr::filter(SAMPLE_ID %in% genelist$gene_id.x) %>%
#   dplyr::filter(CNV==2) %>%
#   dplyr::group_by(SAMPLE_ID) %>%
#   dplyr::select(-sample) %>%
#   dplyr::mutate(homo_amp=sum(CNV)/2) %>%
#   dplyr::select(-CNV) %>%
#   unique() -> common_gene_cnv.homo_amp
# common_gene_cnv %>%
#   dplyr::filter(SAMPLE_ID %in% genelist$gene_id.x) %>%
#   dplyr::filter(CNV==-2) %>%
#   dplyr::group_by(SAMPLE_ID) %>%
#   dplyr::select(-sample) %>%
#   dplyr::mutate(homo_del=sum(CNV)/-2) %>%
#   dplyr::select(-CNV) %>%
#   unique() -> common_gene_cnv.homo_dele
# common_gene_cnv$sample %>% unique() %>% length() -> smaple_n
# 
# common_gene_cnv.homo_amp %>%
#   dplyr::full_join(common_gene_cnv.homo_dele,by="SAMPLE_ID") %>%
#   dplyr::mutate(n=smaple_n) %>%
#   dplyr::mutate(homo_amp=ifelse(!is.na(homo_amp),100*homo_amp/n,0)) %>%
#   dplyr::mutate(homo_del=ifelse(!is.na(homo_del),100*homo_del/n,0)) %>%
#   dplyr::arrange(desc(homo_del)) -> DOWN_commone_targets.CNV_percent
# DOWN_commone_targets.CNV_percent %>%
#   readr::write_tsv(file.path(data_path,"Figure4/Figure5","Down_common_targets.CNV-percent.tsv"))

# draw figure -------------------------------------------------------------

# DOWN_commone_targets.CNV_percent %>%
#   dplyr::filter(homo_del>5) %>%
#   ggplot(aes(x=SAMPLE_ID,y=homo_del,fill=SAMPLE_ID)) +
#   geom_bar(stat = "identity") +
#   xlab("Gene") +
#   ylab("CNV Deletion (%)") +
#   theme(
#     panel.background = element_blank(),
#     panel.border = element_rect(fill='transparent', color='black'),
#     axis.text = element_text(size = 13),
#     axis.title = element_text(size = 15)
#   )
# ggsave(file.path(data_path,"Figure4/Figure5","Down_common_targets.CNV-percent(del_5).pdf"))

luad_cnv %>%
  dplyr::filter(symbol %in% genelist) %>%
  dplyr::select(symbol,a_homo,d_homo) %>%
  tidyr::gather(-symbol,key="type",value="CNV") %>%
  dplyr::mutate(Percent=CNV*100) %>%
  dplyr::mutate(type=ifelse(type=="a_homo","Amplification","Deletion"))-> plot_ready

plot_ready %>%
  readr::write_tsv(file.path("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/","Figure4/Figure5","Down_common_targets.CNV-percent.tsv"))
plot_ready %>%
  dplyr::filter(type=="Deletion") %>%
  # dplyr::filter(type=="Amplification") %>%
  dplyr::filter(CNV > 0.05) -> more_than5_percent

plot_ready %>%
  dplyr::filter(symbol %in% more_than5_percent$symbol) %>%
  ggplot(aes(x=symbol,y=Percent,fill=type)) +
  geom_col(position = "stack",width = 0.5) +
  guides(fill=guide_legend(title = "")) +
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  ylab("Percent (%)") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill='transparent',colour = "black"),
    legend.position = c(0.3,0.9),
    legend.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size = 12)
  ) ->p;p
ggsave(file.path(result_path,"Figure4/Figure5","Down_common_targets.CNV-percent(del_5).pdf"),width = 4,height = 3)
ggsave(file.path(result_path,"Figure4/Figure5","Down_common_targets.CNV-percent(del_5).tiff"),width = 4,height = 3)
ggsave(file.path(result_path,"Figure1","cellcycle.CNV-percent(Ampl_5).pdf"),width = 4,height = 3)
ggsave(file.path(result_path,"Figure1","cellcycle.CNV-percent(Ampl_5).tiff"),width = 4,height = 3)
ggsave(file.path(result_path,"Figure1","ppar.CNV-percent(Ampl_5).pdf"),width = 4,height = 3)
ggsave(file.path(result_path,"Figure1","ppar.CNV-percent(Ampl_5).tiff"),width = 4,height = 3)
