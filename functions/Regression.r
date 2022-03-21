suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/Drug_class_complete/Biogeographic_regions/ppm.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
   make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/drug_plants/edaphic_properties/", 
              help="output file name [default= %default]", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-A' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

 if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-O' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outdir <- opt$out  
  #setwd(wd_location)  
  }

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

library(data.table)
library(openxlsx)
library(tidyr)
library(tibble)
library("ggpubr")
metadata<-read.xlsx(metafile)
countdata<-fread(abundance) 
options(scipen = 999)
countdata <- data.frame(t(countdata[-1]))
header.true <- function(df) {
names(df) <- as.character(unlist(df[1,]))
df[-1,]
}
countdata= header.true(countdata)
countdata$Unmapped=NULL
metadata<-read.xlsx("/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx")
metadata$BARCODE_ID<-paste0("L",metadata$BARCODE_ID)
countdata <- tibble::rownames_to_column(countdata, "BARCODE_ID")
countmeta = merge(metadata,countdata, by = "BARCODE_ID")
sapply(countmeta[30:390], as.numeric)
for (i in colnames(countmeta[15:26])){
for (j in colnames(countmeta[30:390])){
dir.create(paste0("edaphic_properties/",j))
b=cor(as.numeric(countmeta[, i ]), as.numeric(countmeta[, j]), method = c("spearman"), use = "complete.obs")
b= round(b, digits=4)
b= paste(b, "_", sep="")
plotsc=ggscatter(countmeta, x = i, y = j, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          )
ggsave(path=paste(outdir, j, sep="/"), file=paste(b, i, ".pdf", sep=""),plot=plotsc, width = 30, height = 60, limitsize = FALSE, units = "cm", device="pdf")
}}
