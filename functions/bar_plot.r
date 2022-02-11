#stacked bar plot and bar plot.
#option T: excel file containing ppm counts. In this example, the 10 most abundant drug class in italy divided by country. Such file is produced via the DESeq.r function and select_samples.r
#option C: txt file containing the differential abundance analysis results. Such file is produced via the DESeq.r function.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-T", "--tableppm"), type="character", default="/mnt/vol1/projects/LUCAS/functions/ppm_average_terrain_type.xlsx", 
              help="output file name [default= %default]", metavar="character"),
   make_option(c("-C", "--comparison"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/IT_LC_simpl_2018/LC_simpl_2018_Woodland_vs_Grassland.txt",
               help="output file name [default= %default]", metavar="character"),
  make_option(c("-O", "--outputname"), type="character", default="Woodland_vs_Grassland",
	      help="output file name [default= %default]", metavar="character"),
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
if (is.null(opt$outputname)) {
  stop("WARNING: No token specified with '-=' flag.")
} else {  cat ("token is ", opt$outputname, "\n")
  outputname <- opt$outputname  
  }
    library(data.table)
	library(openxlsx)
	library("RColorBrewer")
	library("ggplot2")
 library(tidyr)
#stacked bar plot
full2=read.xlsx(tableppm)
full2 <- full2 %>%
pivot_longer(names_to = "terrain",
values_to = "ppm",
cols = -Drug_Class)
 
graph= ggplot(full2,aes(fill=Drug_Class, y=ppm, x=terrain)) +
geom_bar(position="stack", stat="identity") +
theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))
ggsave("average.pdf",plot=graph, width = 30, units = "cm", device="pdf")

#function to set the first row of data as column names
header.true <- function(df) {
names(df) <- as.character(unlist(df[1,]))
df[-1,]
}

#bar plot
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
ggsave(file=paste( outputname, ".pdf", sep=""),plot=g, width = 30, units = "cm", device="pdf")
