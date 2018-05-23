
# configuration -----------------------------------------------------------

# data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
# de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"
# data_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
data_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
de_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"

# de_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"

# loading data ------------------------------------------------------------

progene.exp <- read.table(file.path(data_path,"NOISeq_DE_TF_cpm_30_FC2_all.exp.xls"),sep = '\t',header = T)

DE_info <- read.table(file.path(data_path,"NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T)

# data manage -------------------------------------------------------------
library(magrittr)
colnames(progene.exp) %>% grep(".01",.) -> tumor.pos
colnames(progene.exp) %>% grep(".11",.) -> normal.pos
colnames(progene.exp)[tumor.pos] <- paste("T",1:length(tumor.pos),sep = "_")
colnames(progene.exp)[normal.pos] <- paste("N",1:length(normal.pos),sep = "_")

progene.exp %>%
  dplyr::as_tibble() %>%
  dplyr::arrange(gene_id) %>%
  tidyr::gather(-gene_id,key="Sample",value="Exp") %>%
  dplyr::arrange(Sample) %>%
  tidyr::spread(key=Sample,value=Exp) %>%
  as.data.frame() ->progene.exp

sample_info <- data.frame(group=substr(colnames(progene.exp)[-1],1,1) %>% as.character()) 
rownames(sample_info) <- colnames(progene.exp)[-1]


DE_info %>%
  dplyr::arrange(gene_id) %>%
  # dplyr::mutate(log2T_mean=log2(case_mean)) %>%
  # dplyr::mutate(log2N_mean=log2(con_mean)) %>%
  dplyr::select(gene_id,log2FC) -> DE_updown_info #log2T_mean,log2N_mean,
rownames(DE_updown_info) <- DE_updown_info$gene_id
gene_info <- as.data.frame(DE_updown_info[,-1],ncol=1)
rownames(gene_info) <- rownames(DE_updown_info)
colnames(gene_info) <- "log2FC"
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
                                                                      # c("#00C5CD","white", "#D15FEE"))
                                      ),
                           width = unit(0.5, "cm"))
draw(gene_anno,1:20)

sample_anno <- HeatmapAnnotation(df = sample_info,
                                 col = list(group=c("T" = "#8C8C8C", "N" = "#FFC1C1")),
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

out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure1"

he = Heatmap(progene.exp.scaled,
             show_row_names = FALSE, 
             show_column_names = FALSE,
             cluster_columns = FALSE,
             top_annotation = sample_anno,
             show_row_dend = FALSE, # whether show row clusters.
             heatmap_legend_param = list(title = c("Experssion")))
pdf(file.path(out_path,"FC2_progene_exp_heatmap.pdf"),width = 6,height = 6)
he+gene_anno
dev.off()
