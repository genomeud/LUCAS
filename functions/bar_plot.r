#stacked bar plot and bar plot

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-T", "--tableppm"), type="character", default="final10IT.xlsx", 
              help="output file name [default= %default]", metavar="character"),
   make_option(c("-C", "--comparison"), type="character", default="LC_simpl_2018_Cropland_vs_Grassland.txt",
               help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USATE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$tableppm)) {
  stop("WARNING: No token specified with '-T' flag.")
} else {  cat ("token is ", opt$tableppm, "\n")
  tableppm <- opt$tableppm  
  }
if (is.null(opt$comparison)) {
  stop("WARNING: No token specified with '-C' flag.")
} else {  cat ("token is ", opt$tableppm, "\n")
  comparison <- opt$comparison  
  }
    library(data.table)
	library(ape)
	library(openxlsx)
	library("RColorBrewer")
	library("pheatmap")
	library("ggplot2")
	library("ggrepel")
 library(tidyr)

full2=read.xlsx(tableppm)
full2 <- full2 %>%
pivot_longer(names_to = "terrain",
values_to = "ppm",
cols = -Drug_Class)
 
graph= ggplot(full2,aes(fill=Drug_Class, y=ppm, x=terrain)) +
geom_bar(position="stack", stat="identity") +
theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))
ggsave("ringo.pdf",plot=graph, width = 30, units = "cm", device="pdf")

browser()

header.true <- function(df) {
names(df) <- as.character(unlist(df[1,]))
df[-1,]
}


full2=read.table(comparison)
full2=header.true(full2)
full2 = subset(full2, padj < 0.05)
full2$stat=full2$padj=full2$lfcSE=full2$pvalue=full2$baseMean=NULL
full2$padj= as.numeric(as.character(full2$log2FoldChange))
g <- ggplot(full2, aes(x = class , y = log2FoldChange)) +
  geom_bar(
    stat = "identity", position = position_stack(),
    color = "white", fill = "lightblue"
  ) +
  coord_flip()
ggsave("averdasadge.pdf",plot=g, width = 30, units = "cm", device="pdf")
