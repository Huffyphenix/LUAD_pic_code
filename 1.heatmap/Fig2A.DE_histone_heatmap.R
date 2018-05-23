# configuration -----------------------------------------------------------
# 
# data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
# de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"
data_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
de_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
gene_info <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TCGA_gene_info"
gene_info <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/TCGA_gene_info"

# data_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
# de_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"


# loading data ------------------------------------------------------------

TF.exp <- read.table(file.path(data_path,"NOISeq_DE_TF_cpm_30_FC2_all.exp.xls"),sep = '\t',header = T)
progene.exp <- read.table(file.path(data_path,"NOISeq_DE_ProGene_FC2_cpm_30.exp.xls"),sep = '\t',header = T)
rbind(TF.exp,progene.exp) -> all_gene_de_exp
library(magrittr)
TF_DE_info <- read.table(file.path(de_path,"NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T) 
progene_DE_info <- read.table(file.path(de_path,"NOISeq_DE_ProGene_FC2_cpm_30"),sep = '\t',header = T) 
rbind(TF_DE_info,progene_DE_info) -> all_DE_info

# histone <- readr::read_tsv("H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/hisone/histone-methylation/all_histone_methylation.idmap",col_names = F)
histone <- c("CBX2","EZH2","CBX7","CBX3","DNMT3A","CBX8","SUV39H2","UHRF1")
tcga_geneinfo <- readr::read_tsv(file.path(gene_info,"TCGA_all_gene_id.txt"))
ncbi_geneinfo_9606 <- readr::read_tsv(file.path(gene_info,"Homo_sapiens.gene_info.filter.20180423download-ncbi"))
ncbi_geneinfo_9606 %>%
  dplyr::filter(Symbol %in% histone) -> histone_geneID_normorlize
tcga_geneinfo %>%
  dplyr::filter(entrez_id %in% histone_geneID_normorlize$GeneID) -> histone_tcga_geneID_normorlize
# data manage -------------------------------------------------------------
all_gene_de_exp %>%
  dplyr::filter(gene_id %in% histone_tcga_geneID_normorlize$gene_id) -> histone_exp
all_DE_info %>%
  dplyr::filter(gene_id %in% histone_tcga_geneID_normorlize$gene_id) -> histone_de_info

library(magrittr)
colnames(histone_exp) %>% grep(".01",.) -> tumor.pos
colnames(histone_exp) %>% grep(".11",.) -> normal.pos
colnames(histone_exp)[tumor.pos] <- paste("T",1:length(tumor.pos),sep = "_")
colnames(histone_exp)[normal.pos] <- paste("N",1:length(normal.pos),sep = "_")

histone_exp %>%
  dplyr::as_tibble() %>%
  tidyr::gather(-gene_id,key="Sample",value="Exp") %>%
  dplyr::arrange(Sample) %>%
  tidyr::spread(key=Sample,value=Exp) %>%
  dplyr::arrange(gene_id) %>%
  as.data.frame() -> histone_exp

sample_info <- data.frame(group=substr(colnames(histone_exp)[-1],1,1) %>% as.character()) 
rownames(sample_info) <- colnames(histone_exp)[-1]


histone_de_info %>%
  # dplyr::mutate(log2T_mean=log2(case_mean)) %>%
  # dplyr::mutate(log2N_mean=log2(con_mean)) %>%
  dplyr::select(gene_id,log2FC) %>%
  dplyr::arrange(log2FC) -> DE_updown_info #log2T_mean,log2N_mean,
rownames(DE_updown_info) <- DE_updown_info$gene_id
gene_info <- as.data.frame(DE_updown_info[,-1],ncol=1)
rownames(gene_info) <- rownames(DE_updown_info)
colnames(gene_info) <- "log2FC"
histone_exp %>%
  dplyr::inner_join(DE_updown_info,by="gene_id") %>%
  dplyr::arrange(log2FC) %>%
  dplyr::select(-log2FC)  -> histone_exp
# draw pic ----------------------------------------------------------------
library(ComplexHeatmap)
# gene row annotation
gene_anno <- rowAnnotation(df=gene_info,
                           col = list(log2FC=circlize::colorRamp2(c(min(gene_info$log2FC),
                                                                    0,
                                                                    max(gene_info$log2FC)),
                                                                  c("green","white","red"))
                                      # log2N_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                      #                                   median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                      #                                   max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                      #                                 c("#00C5CD","white", "#D15FEE")),
                                      # log2T_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                      #                                   median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                      #                                   max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                      #                                 c("#00C5CD","white", "#D15FEE"))
                                      ),
                           annotation_legend_param = list(title = c("log2(FC)")),
                           # annotation_name_side = "bottom",
                           # annotation_name_rot = 180,
                           width = unit(0.5, "cm"))
draw(gene_anno,1:20)

sample_anno <- HeatmapAnnotation(df = sample_info,
                                 col = list(group=c("T" = "#8C8C8C", "N" = c("#FFC1C1"))),
                                 width = unit(0.5, "cm"),
                                 annotation_legend_param = list(title = c("Group")),
                                 name = "Group")
draw(sample_anno,1:118)

# sample column annotation
rownames(histone_exp) <- histone_exp[,1]
histone_exp <- histone_exp[,-1]
histone_exp <- as.matrix(histone_exp)

histone_exp.scaled <- apply(histone_exp,1,scale) %>% t()
# rownames(progene.exp.scaled) <- rownames(progene.exp)
colnames(histone_exp.scaled) <- colnames(histone_exp)

out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2"

he = Heatmap(histone_exp.scaled,
             show_row_names = TRUE, 
             show_column_names = FALSE,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             top_annotation = sample_anno,
             heatmap_legend_param = list(title = c("Scaled Exp")))
pdf(file.path(out_path,"DE_histone_heatmap.pdf"),width = 6,height = 3)
gene_anno+he
dev.off()
