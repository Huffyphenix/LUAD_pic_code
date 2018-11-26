
# configuration -----------------------------------------------------------

# data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/??ͼ"
# de_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/????????"

# data_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/热图"
# de_path <- "F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/差异表达data"
data_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图"
de_path <- "S:/study/ENCODE-TCGA-LUAD/result/热图/20160519.FC2"
# loading data ------------------------------------------------------------

library(magrittr)
progene.exp <- read.table(file.path(data_path,"20160519.FC2","NOISeq_DE_ProGene_FC2_cpm_30.exp.xls"),sep = '\t',header = T)
TF_exp <- read.table(file.path(data_path,"20160519.FC2","NOISeq_DE_TF_cpm_30_FC2_all.exp.xls"),sep = '\t',header = T)
progene.exp %>%
  rbind(TF_exp) -> progene.exp
pro_DE_info <- read.table(file.path(de_path,"NOISeq_DE_ProGene_FC2_cpm_30"),sep = '\t',header = T) %>%
  dplyr::arrange(gene_id)
tf_DE_info <- read.table(file.path(de_path,"NOISeq_DE_TF_FC2_cpm_30"),sep = '\t',header = T) %>%
  dplyr::arrange(gene_id)
pro_DE_info %>%
  rbind(tf_DE_info) %>%
  dplyr::arrange(gene_id) -> DE_info

# data manage -------------------------------------------------------------
library(magrittr)
colnames(progene.exp) %>% grep(".01",.) -> tumor.pos
colnames(progene.exp) %>% grep(".11",.) -> normal.pos
colnames(progene.exp)[tumor.pos] <- paste("T",1:length(tumor.pos),sep = "_")
colnames(progene.exp)[normal.pos] <- paste("N",1:length(normal.pos),sep = "_")

# colnames arrange 
progene.exp %>%
  dplyr::as_tibble() %>%
  tidyr::gather(-gene_id,key="Sample",value="Exp") %>%
  dplyr::arrange(Sample) %>%
  tidyr::spread(key=Sample,value=Exp) %>%
  as.data.frame() ->progene.exp
sample_info <- data.frame(group=substr(colnames(progene.exp)[-1],1,1) %>% as.character()) 
rownames(sample_info) <- colnames(progene.exp)[-1]

progene.exp %>%
  dplyr::as_tibble() %>%
  tidyr::gather(-gene_id,key="Sample",value="Exp") %>%
  dplyr::mutate(Group = substr(Sample,1,1)) %>%
  dplyr::select(-Sample) %>%
  dplyr::group_by(gene_id,Group) %>%
  dplyr::do(
    Mean=sum(.$Exp)/58
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(Mean=scale(Mean)) %>%
  tidyr::spread(Group,Mean) %>%
  dplyr::rename("N_mean"="N","T_mean"="T") %>%
  as.data.frame() -> progene.exp.mean
rownames(progene.exp.mean) <- progene.exp.mean[,1]
progene.exp.mean <- progene.exp.mean[,-1]

DE_info %>%
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
                                                                      c("green","white","red"))),
                                          # log2N_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                          #                               median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                          #                               max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                          #                             c("blue","white","red")),
                                          # log2T_mean=circlize::colorRamp2(c(min(min(gene_info$log2N_mean),min(gene_info$log2T_mean)), 
                                          #                               median(c(gene_info$log2N_mean,gene_info$log2T_mean)),
                                          #                               max(max(gene_info$log2N_mean),max(gene_info$log2T_mean))),
                                          #                             c("blue","white", "red"))
                                          
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

out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/supplymentary"
library(circlize)
he = Heatmap(progene.exp.scaled,
             col = colorRamp2(c(-2, 0, 4), c(c("#00BFFF"), "white", "red")),
        show_row_names = FALSE, 
        show_column_names = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE, # whether show row clusters.
        top_annotation = sample_anno,
        heatmap_legend_param = list(title = c("Expression")))
pdf(file.path(out_path,"FC2_progene_exp_heatmap.pdf"),width = 5,height = 7)
he+gene_anno
dev.off()
