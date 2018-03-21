library(magrittr)
# data path ---------------------------------------------------------------

data_path <- "H:/WD Backup.swstor/MyPC/MDNkNjQ2ZjE0ZTcwNGM0Mz/Volume{3cf9130b-f942-4f48-a322-418d1c20f05f}/study/ENCODE-TCGA-LUAD/result/EZH2分析/甲基化分析/图"

# laod data ---------------------------------------------------------------

cbx2 <- readr::read_tsv(file.path(data_path,"CBX2_promoter.txt"))
ezh2 <- readr::read_tsv(file.path(data_path,"EZH2_promoter.txt"))


# data manage ------------------------------------------------------------
cbx2 %>%
  dplyr::mutate(Promoter=c(1:nrow(cbx2))) %>%
  tidyr::gather(-tag,-Promoter,key="Group",value="Methylation") %>%
  dplyr::mutate(gene="CBX2") -> cbx2_gather
ezh2 %>%
  dplyr::mutate(Promoter=seq(1,nrow(cbx2),10/6)) %>%
  tidyr::gather(-tag,-Promoter,key="Group",value="Methylation") %>%
  dplyr::mutate(gene="EZH2") -> ezh2_gather

rbind(cbx2_gather,ezh2_gather) -> all_methy


# draw pic ----------------------------------------------------------------
library(ggplot2)
all_methy %>%
  ggplot(aes(x=Promoter,y=Methylation,color=Group)) +
  geom_line() +
  facet_grid(~gene) +
  theme(
    panel.background  = element_rect(fill = "white", color = "grey"),
    axis.title.x = element_text()
  )
ggsave("F:/我的坚果云/ENCODE-TCGA-LUAD/Figure/Figure2/Figure2C.Methy_histone.pdf",device = "pdf",width = 4,height = 3)
