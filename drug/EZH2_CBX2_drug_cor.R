library(magrittr)
# data path ---------------------------------------------------------------

drug_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/drug_sensitivity"
gdsc_info <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/DRUG_data/GDSC database"


# load data ---------------------------------------------------------------

EZH2_CBX2_LUAD_IC50_exp <- readr::read_tsv(file.path(drug_path,"EZH2_CBX2_LUAD_IC50_exp.tsv"))
drug_info <- readr::read_tsv(file.path(drug_path,"Screened_Compounds.txt")) %>%
  dplyr::rename("DRUG_ID"="DRUG ID")


# data combination --------------------------------------------------------

EZH2_CBX2_LUAD_IC50_exp %>%
  dplyr::inner_join(drug_info,by="DRUG_ID") -> EZH2_CBX2_LUAD_IC50_exp.druginfo


# calculation -------------------------------------------------------------

fn_cor <- function(data){
  broom::tidy(cor.test(data$exp,data$LN_IC50,method = c("spearman")),
              warning =function(e) 2 ,
              error=function(e) 1) -> tmp.spm
  if(length(tmp.spm)!=1){
    return(tmp.spm)
  }
}
EZH2_CBX2_LUAD_IC50_exp %>%
  tidyr::nest(-ensembl_gene,-DRUG_ID) %>%
  dplyr::group_by(ensembl_gene,DRUG_ID) %>%
  dplyr::mutate(spm=purrr::map(data,fn_cor)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value<=0.05) %>%
  dplyr::left_join(drug_info,by="DRUG_ID") %>%
  dplyr::mutate(symbol=ifelse(ensembl_gene=="ENSG00000173894","CBX2","EZH2")) %>%
  dplyr::select(symbol,`DRUG NAME`,TARGET,`TARGET PATHWAY`,estimate,p.value) -> EZH2_CBX2_cor_drug




# draw pic ----------------------------------------------------------------
EZH2_CBX2_cor_drug %>%
  dplyr::mutate(`TARGET PATHWAY`=ifelse(`TARGET PATHWAY`==TARGET,paste(TARGET,1,sep="_"),`TARGET PATHWAY`)) -> EZH2_CBX2_cor_drug
  

### get text 
EZH2_CBX2_cor_drug %>%
  dplyr::pull(symbol) %>%
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
  dplyr::mutate(drug_target=c("Inhibit","Inhibit","Inhibit","Inhibit","Inhibit","Activate","Inhibit")) %>%
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
out_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure3"
ggsave(file.path(out_path,"EZH2_CBX2_drug_sensiticity.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(out_path,"EZH2_CBX2_drug_sensiticity.tiff"),device = "tiff",height = 6,width = 8)
