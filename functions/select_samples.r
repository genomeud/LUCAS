#extract samples from the ppm table, produced via DESeq2.r, by selecting a variable and a specific value. 
#--abundance : ppm table of sample collected in "italy", produced via DESeq2.r
#--metafile : Metadata of italian samples
#--condition : "LC_simpl_2018", terrain type
#--token : "Woodland", terrain that i want to select
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/IT_LC_simpl_2018/ppm.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/IT_metadata.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="LC_simpl_2018", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-T", "--token"), type="character", default="Woodland", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-A' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

if (is.null(opt$condition)) {
  stop("WARNING: No condition specified with '-C' flag.")
} else {  cat ("condition is ", opt$condition, "\n")
  condition <- opt$condition  
  }

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-M' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

  if (is.null(opt$token)) {
  stop("WARNING: No token specified with '-T' flag.")
} else {  cat ("token is ", opt$token, "\n")
  token <- opt$token  
  }
 
  
      library(data.table)
	library(ape)
	library(openxlsx)
   library("stringr")
  library("dplyr")
  
  metadata<-read.xlsx(metafile)
	countdata<-fread(abundance,data.table=F)
	metadata$BARCODE_ID<-paste0("L",metadata$BARCODE_ID)
 	metadata<-metadata[,c("BARCODE_ID",condition)]
    metadata = lapply(metadata, subset, metadata$LC_simpl_2018 == token)
 	names(metadata)[2]<-"condition"
 	rownames(countdata)<-countdata$Drug_Class
  countdata$Drug_Class = NULL
  #keep only the samples present in metadata
countdata<-countdata[,names(countdata)%in%metadata$BARCODE_ID]
  
write.table(countdata,file=paste( token, ".txt", sep=""), row.names=F, quote=F)
 
