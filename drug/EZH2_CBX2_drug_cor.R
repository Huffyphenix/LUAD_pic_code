library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
# data path in Ezhou---------------------------------------------------------------

drug_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/drug_sensitivity"
gdsc_info <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/DRUG_data/GDSC database"
FFL_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets-180830"
out_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure9"


# data path in Hust -------------------------------------------------------
drug_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/drug_sensitivity"
gdsc_info <- "G:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/DRUG_data/GDSC database"
FFL_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/CBX2_H3K27me3-common-targets/common-targets-180426-new/FFL/EZH2_CBX2_targets-180830"
out_path <- "S:/坚果云/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure9"


# load data ---------------------------------------------------------------

LUAD_IC50_exp <- readr::read_rds(file.path(drug_path,"data_merged","LUAD_cell_line_expression_IC50.rds.gz"))
drug_info <- readr::read_tsv(file.path(drug_path,"Screened_Compounds.txt")) %>%
  dplyr::rename("DRUG_ID"="DRUG ID")

PPAR_TSG <- readr::read_tsv(file.path(FFL_path,"attribute_fc_tsg_multuiple-source.txt")) %>%
  dplyr::filter(!SYMBOL %in% c("APOA1","APOA5"))

genelist <- readr::read_tsv(file.path(FFL_path,"attribute_fc_tsg_multuiple-source.txt")) %>%
  .$SYMBOL %>%
  unique() %>%
  bitr(fromType = "SYMBOL", toType = "ENSEMBL",OrgDb = org.Hs.eg.db)

c("EZH2","CBX2") %>%
  bitr(fromType = "SYMBOL", toType = "ENSEMBL",OrgDb = org.Hs.eg.db) -> EZH2_CBX2
c("E2F1","SOX4") %>%
  bitr(fromType = "SYMBOL", toType = "ENSEMBL",OrgDb = org.Hs.eg.db) -> E2F1_SOX4
# data combination --------------------------------------------------------

genelist %>%
  rbind(EZH2_CBX2) %>%
  rbind(E2F1_SOX4) %>%
  dplyr::rename("ensembl_gene"="ENSEMBL") -> genelist

LUAD_IC50_exp %>%
  dplyr::filter(ensembl_gene %in% c(genelist$ensembl_gene)) %>%
  dplyr::left_join(genelist,by="ensembl_gene") %>%
  tidyr::unnest() %>%
  dplyr::inner_join(drug_info,by="DRUG_ID") -> LUAD_IC50_exp.druginfo


# calculation -------------------------------------------------------------

fn_cor <- function(data){
  broom::tidy(cor.test(data$exp,data$LN_IC50,method = c("spearman")),
              warning =function(e) 2 ,
              error=function(e) 1) -> tmp.spm
  if(length(tmp.spm)!=1){
    return(tmp.spm)
  }
}
LUAD_IC50_exp.druginfo %>%
  tidyr::nest(-SYMBOL,-DRUG_ID) %>%
  dplyr::group_by(SYMBOL,DRUG_ID) %>%
  dplyr::mutate(spm=purrr::map(data,fn_cor)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::inner_join(drug_info,by="DRUG_ID") %>%
  dplyr::select(SYMBOL,`DRUG NAME`,TARGET,`TARGET PATHWAY`,estimate,p.value) -> cor_drug

cor_drug %>%
  readr::write_tsv(file.path(out_path,"drug_correlation_result.txt"))


# draw pic ----------------------------------------------------------------
PPAR_TSG %>% 
  dplyr::filter(Class == "TSG") %>%
  .$SYMBOL -> genelist_to_draw

cor_drug %>%
  dplyr::mutate(`TARGET PATHWAY`=ifelse(`TARGET PATHWAY`==TARGET,paste(TARGET,1,sep="_"),`TARGET PATHWAY`)) %>%
  dplyr::filter(SYMBOL %in% genelist_to_draw)-> EZH2_CBX2_cor_drug
  

###### key correlation 
#### filtered correlation pairs to show
readr::read_tsv(file.path(out_path,"key_correlations.txt"))-> EZH2_CBX2_cor_drug
# line plot ---------------------------------------------------------------


### get text 
EZH2_CBX2_cor_drug %>%
  dplyr::pull(SYMBOL) %>%
  unique() -> gene.text

EZH2_CBX2_cor_drug %>%
  dplyr::pull(`DRUG NAME`) %>%
  unique() -> drug.text

EZH2_CBX2_cor_drug %>%
  dplyr::pull(`TARGET`) %>%
  unique() -> target.text

EZH2_CBX2_cor_drug %>%
  dplyr::pull(`TARGET PATHWAY`) %>%
  unique() -> targetpath.text

c.text <- data.frame(x = 0.5, y = 1, text = "test", type = "test")

for (i in 1:length(drug.text)) {
  data.frame(x = 3, y = i * 2 - 1, text = drug.text[i], type = "drug") -> tmp.text
  print(tmp.text)
  rbind(c.text, tmp.text) -> c.text
}
c.text$y %>% max() / length(target.text) -> t.i
for (i in 1:length(target.text)) {
  data.frame(x = 5, y = t.i * i - 1, text = target.text[i], type = "target") -> tmp.text
  print(tmp.text)
  rbind(c.text, tmp.text) -> c.text
}
c.text$y %>% max() / length(gene.text) -> g.i
for (i in 1:length(gene.text)) {
  data.frame(x = 1, y = i * g.i - 1, text = gene.text[i], type = "gene") -> tmp.text
  rbind(c.text, tmp.text) -> c.text
}

c.text$y %>% max() / length(targetpath.text) -> tp.i
for (i in 1:length(targetpath.text)) {
  data.frame(x = 7, y = i * tp.i - 1, text = targetpath.text[i], type = "targetpath") -> tmp.text
  print(tmp.text)
  rbind(c.text, tmp.text) -> c.text
}


c.text[-1, ] -> c.text

# get seg ----
get_rppa_seg <- function(data,cancer_text) {
  # name <- c("x1","y1","x2","y2","Cancer","Regulation")
  # print(n)
  data[1,1] %>% as.character() -> gene
  data[1,2] %>% as.character() -> drug
  data[1,3] %>% as.character() -> target
  data[1,4] %>% as.character() -> targetpath
  data[1,5] %>% as.numeric() -> diff
  data[1,6] %>% as.character() -> drug_target
  if (diff > 0) {
    line_type <- "Inhibit"
  } else {
    line_type <- "Activate"
  }
  cancer_text %>%
    dplyr::filter(text %in% gene) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(x = x + 0.5) -> g.pos
  
  cancer_text %>%
    dplyr::filter(text %in% drug) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(x = x + 1.5) -> d2.pos
  cancer_text %>%
    dplyr::filter(text %in% drug) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(x = x ) -> d1.pos
  
  cancer_text %>%
    dplyr::filter(text %in% target) %>%
    dplyr::select(x, y)  %>%
    dplyr::mutate(x = x)-> t1.pos
  cancer_text %>%
    dplyr::filter(text %in% target) %>%
    dplyr::select(x, y)  %>%
    dplyr::mutate(x = x + 1.5)-> t2.pos
  
  
  cancer_text %>%
    dplyr::filter(text %in% targetpath) %>%
    dplyr::select(x, y) -> tp.pos
  
  .d_seq_tmp1 <- data.frame(x1 = g.pos$x, y1 = g.pos$y, x2 = d1.pos$x, y2 = d1.pos$y, Pathway = targetpath, Regulation = line_type)
  .d_seq_tmp2 <- data.frame(x1 = d2.pos$x, y1 = d2.pos$y, x2 = t1.pos$x, y2 = t1.pos$y, Pathway = targetpath, Regulation = drug_target)
  .d_seq_tmp3 <- data.frame(x1 = t2.pos$x, y1 = t2.pos$y, x2 = tp.pos$x, y2 = tp.pos$y, Pathway = targetpath, Regulation = "Activate")
  
  rbind(.d_seq_tmp1,.d_seq_tmp2,.d_seq_tmp3) -> .d_seg
  .d_seg$Pathway <- .d_seg$Pathway %>% as.character()
  .d_seg$Regulation <- .d_seg$Regulation %>% as.character()
  tibble::as_tibble(.d_seg)
}
EZH2_CBX2_cor_drug %>%
  dplyr::select(-p.value) %>%
  dplyr::mutate(drug_target=c("Activate")) %>%
  dplyr::mutate(n=1:nrow(EZH2_CBX2_cor_drug)) %>%
  tidyr::nest(-n) %>%
  dplyr::group_by(n) %>%
  dplyr::mutate(seg=purrr::map(data,cancer_text=c.text,.f=get_rppa_seg)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n,-data) %>%
  tidyr::unnest() ->plot_seg

EZH2_CBX2_cor_drug %>%
  dplyr::select(-estimate,-p.value) %>%
  dplyr::mutate(Pathway=`TARGET PATHWAY`) %>%
  tidyr::gather(-Pathway,key="key",value="text") %>%
  dplyr::select(-key) %>%
  dplyr::inner_join(c.text,by="text") -> plotready.text
  
  
# draw pic ----
library(ggplot2)
ggplot() -> p
for (pathway in plot_seg$Pathway %>% unique()) {
  # cancers="LUSC"
  plot_seg %>%
    dplyr::filter(Pathway == pathway) -> data
  curvature <- runif(1, 0.1, 0.3)
  p +
    geom_segment(
      data = data, mapping = aes(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2,
        colour = Pathway,
        linetype = Regulation
        
      ),
      size=1
      # colour = "red",
      # curvature = curvature,
      # arrow = arrow(length = unit(0.03, "npc"))
    ) -> p
}
plotready.text %>%
  dplyr::filter(type == "drug") -> drug.text
plotready.text %>%
  dplyr::filter(type == "target") -> target.text
plotready.text %>%
  dplyr::filter(type == "gene") -> gene.text
plotready.text %>%
  dplyr::filter(type == "targetpath") -> targetpath.text
p +
  guides(color = FALSE) +
  geom_text(
    data = drug.text,
    mapping = aes(x = x, y = y, label = text, color = Pathway),
    hjust = 0,
    size = 4
  ) +
  geom_text(
    data = target.text,
    mapping = aes(x = x , y = y, label = text, color = Pathway),
    hjust = 0,
    size = 4
  ) +
  geom_text(
    data = gene.text,
    mapping = aes(x = x, y = y, label = text),
    hjust = 0,
    size = 4
  ) +
  geom_text(
    data = targetpath.text,
    mapping = aes(x = x, y = y, label = text, color = Pathway),
    hjust = 0,
    size = 4
  ) +
  expand_limits(x = c(2, 10)) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # text = element_text(size=5),
    plot.title = element_text(hjust = 0.5, size = 25),
    plot.margin = rep(unit(0, "null"), 4),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.25, "cm"),
    legend.title = element_text(size = 12),
    text = element_text(size = 20)
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Drug sensitiviy") -> p2;p2
ggsave(file.path(out_path,"EZH2_CBX2_drug_sensiticity.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(out_path,"EZH2_CBX2_drug_sensiticity.tiff"),device = "tiff",height = 6,width = 8)

ggsave(file.path(out_path,"E2F1_SOX4_drug_sensiticity.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(out_path,"E2F1_SOX4_drug_sensiticity.tiff"),device = "tiff",height = 6,width = 8)

ggsave(file.path(out_path,"keys_drug_sensiticity.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(out_path,"keys_drug_sensiticity.tiff"),device = "tiff",height = 6,width = 8)
# point plot --------------------------------------------------------------
EZH2_CBX2_cor_drug %>%
  dplyr::filter(abs(estimate)>=0.3) -> EZH2_CBX2_cor_drug

EZH2_CBX2_cor_drug %>%
  dplyr::arrange(`TARGET PATHWAY`,estimate) %>%
  dplyr::select(`DRUG NAME`,`TARGET PATHWAY`) %>%
  unique() -> drug_rank
EZH2_CBX2_cor_drug %>%
  dplyr::group_by(`TARGET PATHWAY`) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(`TARGET PATHWAY`,n) %>%
  dplyr::arrange(`TARGET PATHWAY`) %>%
  unique() %>%
  dplyr::ungroup() -> n_color

y_coordinate <- vector(length = nrow(n_color)) %>% as.numeric()
for(i in 1:nrow(n_color)){
  if(i==1){
    y_coordinate[i] <- n_color$n[i] %>% as.numeric()
  }else{
    y_coordinate[i] <- y_coordinate[i-1]+n_color$n[i]
  }
}
n_color %>%
  dplyr::mutate(color = ggthemes::gdocs_pal()(nrow(n_color))) %>%
  dplyr::mutate(n_r=1:nrow(n_color)) %>%
  dplyr::mutate(y=(y_coordinate-n/2),x=4) -> pathway_color ## change x when num of gene set changes
drug_rank %>%
  dplyr::left_join(pathway_color,by="TARGET PATHWAY") -> drug_rank
drug_rank %>%
  dplyr::select(-`DRUG NAME`) %>%
  dplyr::distinct() -> text_added
EZH2_CBX2_cor_drug %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::mutate(sum=sum(estimate)) %>%
  dplyr::ungroup() %>%
  dplyr::select(SYMBOL,sum) %>%
  unique() %>%
  dplyr::arrange(sum) -> gene_rank

library(ggplot2)
EZH2_CBX2_cor_drug %>%
  dplyr::mutate(log10p= -log10(p.value)) %>%
  ggplot(aes(x=SYMBOL,y=`DRUG NAME`)) +
  geom_point(aes(size=log10p,color=estimate)) +
  scale_x_discrete(limits = gene_rank$SYMBOL, expand = c(0.012,0.012)) +
  scale_y_discrete(limits = drug_rank$`DRUG NAME`, expand = c(0.012,0.012), position = "right") +
  scale_color_gradient2(
    name = "Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "-log10(p value)"
  ) +
  geom_text(data = text_added,
            aes(x=x,y=y,label=`TARGET PATHWAY`),color=text_added$color) +
  theme(
    panel.background = element_rect(color = "black", fill = "white", size = 0.1),
    panel.grid=element_line(colour="grey",linetype="dashed"),
    panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
    
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10, color = drug_rank$color),
    
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    # legend.key.width = unit(1,"cm"),
    # legend.key.heigh = unit(0.3,"cm"),
    legend.key = element_rect(fill="white",colour = "black")
  ) + guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  ) -> p;p
ggsave(file.path(out_path,"PPAR_drug_sensiticity.pdf"),device = "pdf",height = 5,width = 5)
ggsave(file.path(out_path,"PPAR_drug_sensiticity.tiff"),device = "tiff",height = 5,width = 5)
ggsave(file.path(out_path,"TSG_drug_sensiticity.pdf"),device = "pdf",height = 10,width = 5)
ggsave(file.path(out_path,"TSG_drug_sensiticity.tiff"),device = "tiff",height = 10,width = 5)
