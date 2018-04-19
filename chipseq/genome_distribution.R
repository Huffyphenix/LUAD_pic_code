.libPaths("F:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{7a27e707-64db-4391-94fd-a8b51e3df0b4}/software/R/R-3.4.1/library")

# source("http://bioconductor.org/biocLite.R")
# biocLite("ChIPseeker")
# biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)

# data path ---------------------------------------------------------------

data_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets"
# load data ---------------------------------------------------------------


EZH2_treat <- readPeakFile(file.path(data_path,"H3K27me3_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-treat.2"))
EZH2_control <- readPeakFile(file.path(data_path,"H3K27me3_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-control.2"))
CBX2_treat <-  readPeakFile(file.path(data_path,"CBX2_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-treat.2"))
CBX2_control<-  readPeakFile(file.path(data_path,"CBX2_pva_3-peaks_OVERLAP_GRCh38-hg38_TSS-5kb-control.2"))

covplot(EZH2_treat, weightCol="pileup",lower = 0.001)
covplot(EZH2_control, weightCol="pileup",lower = 0.001)
covplot(CBX2_treat, weightCol="pileup",lower = .001)
covplot(CBX2_control, weightCol="pileup",lower = .001)

# Profile of ChIP peaks binding to TSS regions ----------------------------
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)

tagMatrix_ET <- getTagMatrix(EZH2_treat, windows=promoter)
tagMatrix_EC <- getTagMatrix(EZH2_control, windows=promoter)
tagMatrix_CC <- getTagMatrix(CBX2_control, windows=promoter)
tagMatrix_CT <- getTagMatrix(CBX2_control, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-5000, 5000), color="red")
peakHeatmap(EZH2_treat, TxDb=txdb, upstream=5000, downstream=5000, color="red")
plotAvgProf(tagMatrix_ET, xlim=c(-5000, 5000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix_EC, xlim=c(-5000, 5000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
