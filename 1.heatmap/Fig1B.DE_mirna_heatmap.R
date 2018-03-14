
# configuration -----------------------------------------------------------

data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达"
# loading data ------------------------------------------------------------

progene.exp <- read.table(file.path(data_path,"NOISeq_DE_mirna_FC2_cpm30.exp.mirnaid.xls"),sep = '\t',header = T)

DE_info <- read.table(file.path(de_path,"FC2","NOISeq_DE_mirna_FC2_cpm30.mirnaid"),sep = '\t',header = T)

# data manage -------------------------------------------------------------
library(magrittr)
colnames(progene.exp) %>% grep(".01",.) -> tumor.pos
colnames(progene.exp) %>% grep(".11",.) -> normal.pos
colnames(progene.exp)[tumor.pos] <- paste("T",1:length(tumor.pos),sep = "_")
colnames(progene.exp)[normal.pos] <- paste("N",1:length(normal.pos),sep = "_")
sample_info <- data.frame(group=substr(colnames(progene.exp)[-1],1,1) %>% as.character()) 
rownames(sample_info) <- colnames(progene.exp)[-1]


DE_info %>%
  dplyr::mutate(log2T_mean=log2(case_mean)) %>%
  dplyr::mutate(log2N_mean=log2(con_mean)) %>%
  dplyr::select(mirna_id,log2T_mean,log2N_mean,log2FC) -> DE_updown_info
rownames(DE_updown_info) <- DE_updown_info$mirna_id
gene_info <- DE_updown_info[,-1]
# draw pic ----------------------------------------------------------------

library(ComplexHeatmap)
# gene row annotation
gene_anno <- rowAnnotation(df=gene_info,
                           col = list(log2FC=circlize::colorRamp2(c(min(gene_info$log2FC),
                                                                    0,
                                                                    max(gene_info$log2FC)),
                                                                  c("green","white","red")),
                                      log2N_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                                                        median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                                                        max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                                                      c("blue","white","red")),
                                      log2T_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                                                        median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                                                        max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                                                      c("blue","white", "red"))),
                           width = unit(1.5, "cm"),
                           gap = unit(c(1), "mm"))
draw(gene_anno,1:20)

sample_anno <- HeatmapAnnotation(df = sample_info,
                                 col = list(group=c("T" = "#D1EEEE", "N" = "#FFEFDB")),
                                 width = unit(0.5, "cm"),
                                 name = "Group")
draw(sample_anno,1:118)
# sample column annotation
rownames(progene.exp) <- progene.exp[,1]
progene.exp <- progene.exp[,-1]
progene.exp <- as.matrix(progene.exp)

progene.exp.scaled <- apply(progene.exp,1,scale) %>% t()
# rownames(progene.exp.scaled) <- rownames(progene.exp)
colnames(progene.exp.scaled) <- colnames(progene.exp)

out_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure1"

he = Heatmap(progene.exp.scaled,
             show_row_names = FALSE, 
             show_column_names = FALSE,
             cluster_columns = TRUE,
             top_annotation = sample_anno,
             heatmap_legend_param = list(title = c("Experssion")))
pdf(file.path(out_path,"Figure1B.FC2_mirna_exp_heatmap.pdf"),width = 6,height = 4)
he+gene_anno
dev.off()
