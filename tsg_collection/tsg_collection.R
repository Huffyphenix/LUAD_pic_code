data_path <- "F:/我的坚果云/ENCODE-TCGA-LUAD/TS and oncogene source"


# load data ----------------------------------------------------------------

bushman_onco <- readr::read_tsv(file.path(data_path,"Bushman Lab","allOnco_Feb2017.tsv")) %>%
  dplyr::mutate(source="Brushman") %>%
  dplyr::select(symbol,geneID,source) 
oncogene_database <- readr::read_tsv(file.path(data_path,"oncogene database","oncogene_database-all_the_human_oncogenes.txt")) %>%
  dplyr::mutate(source="oncogene_database") %>%
  dplyr::rename("symbol"="OncogeneName","geneID"="OncogeneID") %>%
  dplyr::select(symbol,geneID,source) 
sanger_census_onco <- readr::read_tsv(file.path(data_path,"Sanger Cancer gene cencus","Tier1_Census_allWed Apr 11 03_17_29 2018.tsv")) %>%
  dplyr::filter(`Role in Cancer` %in% c("oncogene","oncogene, fusion")) %>%
  dplyr::mutate(source="sanger_census") %>%
  dplyr::rename("symbol"="Gene Symbol","geneID"="Entrez GeneId") %>%
  dplyr::select(symbol,geneID,source) 
sanger_census_TSG <- readr::read_tsv(file.path(data_path,"Sanger Cancer gene cencus","Tier1_Census_allWed Apr 11 03_17_29 2018.tsv")) %>%
  dplyr::filter(`Role in Cancer` %in% c("TSG","TSG, fusion")) %>%
  dplyr::mutate(source="sanger_census") %>%
  dplyr::rename("symbol"="Gene Symbol","geneID"="Entrez GeneId") %>%
  dplyr::select(symbol,geneID,source) 
sanger_census_twoside <- readr::read_tsv(file.path(data_path,"Sanger Cancer gene cencus","Tier1_Census_allWed Apr 11 03_17_29 2018.tsv")) %>%
  dplyr::filter(`Role in Cancer` %in% c("oncogene, TSG","oncogene, TSG, fusion")) %>%
  dplyr::mutate(source="sanger_census") %>%
  dplyr::rename("symbol"="Gene Symbol","geneID"="Entrez GeneId") %>%
  dplyr::select(symbol,geneID,source) 

Vogelstein_onco <- readr::read_tsv(file.path(data_path,"Table of oncogenes and tumor suppressor genes from Vogelstein et al. 2013","Driver genes affected by subtle mutations.txt")) %>%
  dplyr::filter(`Classification*` %in% "Oncogene") %>%
  dplyr::mutate(source="Vogelstein.2013") %>%
  dplyr::rename("symbol"="Gene Symbol") %>%
  dplyr::select(symbol,source) 
Vogelstein_onco$symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>%
  dplyr::rename("symbol"="SYMBOL","geneID"="ENTREZID") %>%
  dplyr::inner_join(Vogelstein_onco,by="symbol") -> Vogelstein_onco

Vogelstein_TSG <- readr::read_tsv(file.path(data_path,"Table of oncogenes and tumor suppressor genes from Vogelstein et al. 2013","Driver genes affected by subtle mutations.txt")) %>%
  dplyr::filter(`Classification*` %in% "TSG")  %>%
  dplyr::mutate(source="Vogelstein.2013") %>%
  dplyr::rename("symbol"="Gene Symbol") %>%
  dplyr::select(symbol,source) 
Vogelstein_TSG$symbol %>%
  bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>%
  dplyr::rename("symbol"="SYMBOL","geneID"="ENTREZID") %>%
  dplyr::inner_join(Vogelstein_TSG,by="symbol") -> Vogelstein_TSG
TSGene <- readr::read_tsv(file.path(data_path,"TSGene","LUAD_downregulate_TSG_in_TCGA-TSGeneDatabase.txt"))  %>%
  dplyr::mutate(source="TSGene") %>%
  dplyr::rename("symbol"="GeneSymbol","geneID"="GeneID") %>%
  dplyr::select(symbol,geneID,source) 
UniProt_onco <- readr::read_tsv(file.path(data_path,"UniProt","uniprot-_proto+oncogene_-filtered-organism%3A_Homo+sapiens+%28Human%29+%5B--.tab")) %>%
  dplyr::mutate(source="UniProt") %>%
  tidyr::separate(`Entry name`,c("symbol","species"),"_") %>%
  dplyr::select(Entry,symbol,source) 
UniProt_TSG <- readr::read_tsv(file.path(data_path,"UniProt","uniprot-tumor+suppressor-filtered-organism%3A_Homo+sapiens+%28Human%29+%5B9606--.tab")) %>%
  dplyr::mutate(source="UniProt") %>%
  tidyr::separate(`Entry name`,c("symbol","species"),"_") %>%
  dplyr::select(Entry,symbol,source) 

library(clusterProfiler)
library(org.Hs.eg.db)
UniProt_onco$Entry %>%
  bitr(fromType = "UNIPROT",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>%
  dplyr::rename("Entry"="UNIPROT","geneID"="ENTREZID") %>%
  dplyr::inner_join(UniProt_onco,by="Entry")%>%
  dplyr::select(-Entry,-symbol) ->  UniProt_onco.idfit
UniProt_TSG$Entry %>%
  bitr(fromType = "UNIPROT",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>%
  dplyr::rename("Entry"="UNIPROT","geneID"="ENTREZID") %>%
  dplyr::inner_join(UniProt_TSG,by="Entry")%>%
  dplyr::select(-Entry,-symbol) ->  UniProt_TSG.idfit


# gene combine ------------------------------------------------------------

## oncogene --------------------------------------------------------------
bushman_onco %>%
  rbind(oncogene_database) %>%
  rbind(sanger_census_onco) %>%
  rbind(Vogelstein_onco) %>%
  dplyr::select(-symbol) %>%
  rbind(UniProt_onco.idfit) %>%
  unique() -> all_oncogene
## TS Gene --------------------------------------------------------------
sanger_census_TSG %>%
  rbind(Vogelstein_TSG) %>%
  rbind(TSGene) %>%
  dplyr::select(-symbol) %>%
  rbind(UniProt_TSG.idfit) %>%
  unique() -> all_TSG


### fileter ---------------------------------------------------------
# all_oncogene$geneID %>% table() %>% .[.>=2] %>% names() -> oncogene_at_least_2_source
# all_TSG$geneID %>% table() %>% .[.>=2] %>% names() -> TSG_at_least_2_source
# 
# oncogene_at_least_2_source %>%
#   bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db") -> oncogene_at_least_2_source
# TSG_at_least_2_source %>%
#   bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db") -> TSG_at_least_2_source

# ### confused gene filter --------------------------------------------------------------
# oncogene_at_least_2_source %>%
#   dplyr::inner_join(TSG_at_least_2_source,by="ENTREZID") -> confused_genes
all_oncogene %>%
  dplyr::inner_join(all_TSG,by="geneID") %>%
  dplyr::select(geneID) %>%
  unique() %>%
  .$geneID %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db") -> confused_genes


# ID convertion -----------------------------------------------------------

all_oncogene %>%
  .$geneID %>%
  unique() %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db") -> all_oncogene.symbol

all_TSG %>%
  .$geneID %>%
  unique() %>%
  bitr(fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db") -> all_TSG.symbol


# ###source clear ---------------------------------------------------------

library(plyr)
all_oncogene %>%
  dplyr::rename("ENTREZID"="geneID") %>%
  dplyr::inner_join(all_oncogene.symbol,by="ENTREZID") %>%
  unique() %>%
  ddply(.(ENTREZID,SYMBOL), summarize,
        source=paste(source,collapse=",")) -> oncogene.source_clear

all_TSG %>%
  dplyr::rename("ENTREZID"="geneID") %>%
  dplyr::inner_join(all_TSG.symbol,by="ENTREZID") %>%
  unique() %>%
  ddply(.(ENTREZID,SYMBOL), summarize,
        source=paste(source,collapse=",")) -> TSG.source_clear

# #### confused gene filter------------------------------------------------------

oncogene.source_clear %>%
  dplyr::inner_join(TSG.source_clear,by="ENTREZID") %>%
  dplyr::rename("symbol"="SYMBOL.x","Treat as oncogene"="source.x","Treat as TSG"="source.y") %>%
  dplyr::select(-SYMBOL.y) -> confused_genes.source_clear
sanger_census_twoside %>%
  dplyr::rename("ENTREZID"="geneID") %>%
  dplyr::mutate(`Treat as oncogene` ="sanger_census") %>%
  dplyr::mutate(`Treat as TSG` ="sanger_census") %>%
  dplyr::select(-source) %>%
  rbind(confused_genes.source_clear) %>%
  plyr::ddply(.(symbol,ENTREZID),summarize,
              `Treat as oncogene`=paste(`Treat as oncogene`,collapse=","),
              `Treat as TSG`=paste(`Treat as TSG`,collapse=",")) %>%
  dplyr::filter(! `Treat as oncogene`=="sanger_census")

confused_genes.source_clear %>%
  readr::write_tsv(file.path(data_path,"confuse_gene.source_clear(at least one evidence).tsv"))
TSG.source_clear %>%
  dplyr::filter(!ENTREZID %in% confused_genes.source_clear$ENTREZID) %>%
  readr::write_tsv(file.path(data_path,"TSG.source_clear(at least one evidence-no confuse).tsv"))
oncogene.source_clear %>%
  dplyr::filter(!ENTREZID %in% confused_genes.source_clear$ENTREZID) %>%
  readr::write_tsv(file.path(data_path,"oncogene.source_clear(at least one evidence-no confuse).tsv"))
