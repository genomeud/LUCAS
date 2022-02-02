	# Run with --help flag for help.
# Modified 12/30/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/IT_LC_simpl_2018/ppm.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/IT_metadata.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-R", "--readfile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/readpingas.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/ES_LC_simpl_2018/", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="LC_simpl_2018", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-D", "--value"), type="character", default="grassland", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-G", "--graphdir"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/ES_LC_simpl_2018/", 
              help="output file name [default= %default]", metavar="character")
  # make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              # help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

if (is.null(opt$condition)) {
  stop("WARNING: No condition specified with '-C' flag.")
} else {  cat ("condition is ", opt$condition, "\n")
  condition <- opt$condition  
  }

if (is.null(opt$metafile)) {
  stop("WARNING: No metafile specified with '-V' flag.")
} else {  cat ("metafile is ", opt$metafile, "\n")
  metafile <- opt$metafile  
  }

if (is.null(opt$readfile)) {
  stop("WARNING: No readfile specified with '-V' flag.")
} else {  cat ("readfile is ", opt$readfile, "\n")
  readfile <- opt$readfile  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outdir <- opt$out  
  #setwd(wd_location)  
  }

  if (is.null(opt$graphdir)) {
  stop("WARNING: No graphdir specified with '-G' flag.")
} else {  cat ("graphdir is ", opt$graphdir, "\n")
  graphdir <- opt$graphdir  
  }
 
   if (is.null(opt$graphdir)) {
  stop("WARNING: No graphdir specified with '-D' flag.")
} else {  cat ("graphdir is ", opt$graphdir, "\n")
  graphdir <- opt$graphdir  
  } 
  
     library(DESeq2)
    library(data.table)
	library(ape)
	library(openxlsx)
	library("RColorBrewer")
	library("pheatmap")
	library("ggplot2")
	library("ggrepel")
   library("stringr")
  library("dplyr")
  
  metadata<-read.xlsx(metafile)
	countdata<-fread(abundance,data.table=F)
	metadata$BARCODE_ID<-paste0("L",metadata$BARCODE_ID)
	#Condition based on soil type: LC_simpl_2018
	bigmeta<-metadata
 	metadata<-metadata[,c("BARCODE_ID",condition)]
  metadata = lapply(metadata, subset, metadata$LC_simpl_2018 == "Woodland")
    
 	names(metadata)[2]<-"condition"
	#countdata$`ARO|Name|Drug_Class|Antibiotic|Resistance_Mechanism` = str_c(countdata$ARO, '|' , countdata$Name, '|' , countdata$Drug_Class, '|' , countdata$Antibiotic, '|' , countdata$Resistance_Mechanism)
  #countdata=countdata %>% relocate(`ARO|Name|Drug_Class|Antibiotic|Resistance_Mechanism`)
	#rownames(countdata)<-countdata$`ARO|Name|Drug_Class|Antibiotic|Resistance_Mechanism`
	#countdata$ARO<-countdata$Name<-countdata$Drug_Class<-countdata$Antibiotic<-countdata$Resistance_Mechanism<-NULL
 	rownames(countdata)<-countdata$Drug_Class
  countdata$Drug_Class = NULL
  #keep only the samples present in metadata
 	countdata<-countdata[,names(countdata)%in%metadata$BARCODE_ID]
write.table(countdata,"Woodland", quote=F)
