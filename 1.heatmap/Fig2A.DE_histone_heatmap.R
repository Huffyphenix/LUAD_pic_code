# configuration -----------------------------------------------------------

data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"


# loading data ------------------------------------------------------------

TF.exp <- read.table(file.path(data_path,"NOISeq_DE_TF_cpm_30_FC2_all.exp.xls"),sep = '\t',header = T)
progene.exp <- read.table(file.path(data_path,"NOISeq_DE_ProGene_FC2_cpm_30.exp.xls"),sep = '\t',header = T)
rbind(TF.exp,progene.exp) -> all_gene_de_exp

TF_DE_info <- read.table(file.path(de_path,"FC2","NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T) %>%
  dplyr::rename("Gene_id"="gene_id")
progene_DE_info <- read.table(file.path(de_path,"FC2","NOISeq_DE_ProGene_FC2_cpm_30"),sep = '\t',header = T) 
rbind(TF_DE_info,progene_DE_info) -> all_DE_info

histone <- readr::read_tsv("H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/hisone/histone-methylation/all_histone_methylation.idmap.symbol") %>%
  t() %>% as.character()

# data manage -------------------------------------------------------------
all_gene_de_exp %>%
  dplyr::filter(gene_id %in% histone) -> histone_exp
all_DE_info %>%
  dplyr::filter(Gene_id %in% histone) -> histone_de_info

library(magrittr)
colnames(histone_exp) %>% grep(".01",.) -> tumor.pos
colnames(histone_exp) %>% grep(".11",.) -> normal.pos
colnames(histone_exp)[tumor.pos] <- paste("T",1:length(tumor.pos),sep = "_")
colnames(histone_exp)[normal.pos] <- paste("N",1:length(normal.pos),sep = "_")
sample_info <- data.frame(group=substr(colnames(histone_exp)[-1],1,1) %>% as.character()) 
rownames(sample_info) <- colnames(histone_exp)[-1]


histone_de_info %>%
  dplyr::mutate(log2T_mean=log2(case_mean)) %>%
  dplyr::mutate(log2N_mean=log2(con_mean)) %>%
  dplyr::select(Gene_id,log2T_mean,log2N_mean,log2FC) -> DE_updown_info
rownames(DE_updown_info) <- DE_updown_info$Gene_id
gene_info <- DE_updown_info[,-1]

# draw pic ----------------------------------------------------------------
library(ComplexHeatmap)
# gene row annotation
gene_anno <- rowAnnotation(df=gene_info,
                           col = list(log2FC=circlize::colorRamp2(c(min(gene_info$log2FC),
                                                                    0,
                                                                    max(gene_info$log2FC)),
                                                                  c("#6495ED","white","#FF7F24")),
                                      log2N_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                                                        median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                                                        max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                                                      c("#00C5CD","white", "#D15FEE")),
                                      log2T_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                                                        median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                                                        max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                                                      c("#00C5CD","white", "#D15FEE"))),
                           width = unit(1.5, "cm"),
                           gap = unit(c(1), "mm"))
draw(gene_anno,1:20)

sample_anno <- HeatmapAnnotation(df = sample_info,
                                 col = list(group=c("T" = "#8C8C8C", "N" = "#FFFAFA")),
                                 width = unit(0.5, "cm"),
                                 name = "Group")
draw(sample_anno,1:118)

# sample column annotation
rownames(histone_exp) <- histone_exp[,1]
histone_exp <- histone_exp[,-1]
histone_exp <- as.matrix(histone_exp)

histone_exp.scaled <- apply(histone_exp,1,scale) %>% t()
# rownames(progene.exp.scaled) <- rownames(progene.exp)
colnames(histone_exp.scaled) <- colnames(histone_exp)

out_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2"

he = Heatmap(histone_exp.scaled,
             show_row_names = TRUE, 
             show_column_names = FALSE,
             cluster_columns = TRUE,
             top_annotation = sample_anno,
             heatmap_legend_param = list(title = c("Experssion")))
pdf(file.path(out_path,"DE_histone_heatmap.pdf"),width = 5,height = 5)
he+gene_anno
dev.off()
