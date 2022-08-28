suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/Drug_class_complete/Biogeographic_regions/ppm.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="Country", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/stacked_Drug_class/", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-S", "--second_condition"), type="character", default="LC_simpl_2018", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-T", "--token"), type="character", default="Grassland", 
              help="output file name [default= %default]", metavar="character")
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


if (is.null(opt$condition)) {
  stop("WARNING: No condition specified with '-C' flag.")
} else {  cat ("condition is ", opt$condition, "\n")
  condition <- opt$condition  
  }
  
if (is.null(opt$second_condition)) {
  stop("WARNING: No second_condition specified with '-S' flag.")
} else {  cat ("second_condition is ", opt$second_condition, "\n")
  second_condition <- opt$second_condition  
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
  library(tidyr)
  library(tibble)
  metadata<-read.xlsx(metafile)
	countdata<-fread(abundance,data.table=F)
	metadata$BARCODE_ID<-paste0("L",metadata$BARCODE_ID)
  uni <- distinct(metadata, Country)
  uni <- uni %>% drop_na(Country)
  dir.create(paste0(outdir,condition))
  mydir<-paste0(outdir,condition,"/")
  for (i in 1:length(rownames(uni)))
  {
  paese = uni[i, condition]
  dir.create(paste0(mydir,paese))
  sdir=paste0(mydir,paese,"/")
  metacycle <- metadata[metadata[ , condition] == paese, ]
  metacycle<-metacycle[,c("BARCODE_ID",second_condition)]
  names(metacycle)[2]<- second_condition
  uniz <- distinct(metacycle, LC_simpl_2018)
  uniz <- uniz %>% drop_na(LC_simpl_2018)
  print(paese)
  IDs=fread("/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/Drug_class_complete/LC_simpl_2018/ppm.txt", data.table=F)
  IDs[, 2:628]=NULL
  row_remove <- c("Unmapped")
  unmapped_removed <- IDs[!(row.names(IDs) %in% row_remove),]
  for (j in 1:length(rownames(uniz))) 
  {
  countcycle <- countdata
  suolo = uniz[j, second_condition]
  print(suolo)
  intern_cycle = lapply(metacycle, subset, metacycle[ , second_condition] == suolo)
    #rownames(countdata)<-countdata$`ARO|Name|Drug_Class|Antibiotic|Resistance_Mechanism`
  #countdata$`ARO|Name|Drug_Class|Antibiotic|Resistance_Mechanism` = NULL
  #rownames(countcycle)<-countcycle$Name
  #countcycle$Name=NULL
  rownames(countcycle)<-countcycle$Drug_Class
  countcycle$Drug_Class=NULL
    #keep only the samples present in metadata
  countcycle<-countcycle[,names(countcycle)%in%intern_cycle$BARCODE_ID]
  
      #if (length(colnames(countcycle) > 3) == TRUE)
  long = length(colnames(countcycle))
  if (long > 1)
  {
  row_remove <- c("Unmapped")
  unmapped_removed <- countcycle[!(row.names(countcycle) %in% row_remove),]
  rsfile<-paste(sdir,suolo,".txt",sep="")
  write.table(unmapped_removed,rsfile, quote=F)
  countcycle[ , suolo] <- rowMeans(countcycle, na.rm=TRUE)
  countcycle <- select(countcycle, -contains("L", ignore.case = FALSE))
  countcycle <- tibble::rownames_to_column(countcycle, "class")
  IDs <- merge (IDs,countcycle, by = "class")
  rsfile<-paste(sdir,suolo,"_average.txt",sep="")
  write.table(countcycle,rsfile, quote=F, row.names = F)
  
  }
 
}
print("graph")
IDs_length = length(colnames(IDs))
if (IDs_length > 1)
{
IDs <- IDs %>% mutate(total=rowSums(select_if(., is.numeric)))
  IDs <- IDs[order(-IDs$total),]
  IDs <- head(IDs,10)
  most <- IDs
  most <- most %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
  most <- most[order(-most$total),]
  most <- head(most,1)
  print(most)
  testit <- function(x)
  {
    p1 <- proc.time()
    Sys.sleep(x)
    proc.time() - p1 # The cpu usage should be negligible
  }
  testit(1)
  IDs$total <- NULL
  IDs <- IDs %>%
  pivot_longer(names_to = "Soil_type",
  values_to = "ppm",
  cols = -class)
  graph= ggplot(IDs,aes(fill=class, y=ppm, x=Soil_type)) +
  ggtitle(paese) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))
  resfile<-paste(sdir,paese,"_stacked.pdf",sep="")
  ggsave(resfile,plot=graph, width = 40, units = "cm", device="pdf")
  
  #write.table(countcycle,resfile, quote=F)
  }
}


 
