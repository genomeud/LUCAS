# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/grouped_table/raw_plants_drug/average/average_raw_drug_plants.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-R", "--readfile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/countreads.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/grouped_table/comparison", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="LC_simpl_2018", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-G", "--graphdir"), type="character", default="/mnt/vol1/projects/LUCAS/grouped_table/raw_plants_drug/average/", 
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

runDEseq<-function() 
{
    library(DESeq2)
    library(data.table)
	library(ape)
	library(openxlsx)
	library("RColorBrewer")
	library("pheatmap")
	library("ggplot2")
	library("ggrepel")
	#Create a prefix for the output graph. 
	#WARNING: assumes that the only "." in the filename is the one before file extension and gets as name everything before the "."
	prefname<-unlist(lapply(strsplit(basename(abundance),"\\."),"[",1))
	countdata<-fread(abundance,data.table=F)
	#The first column brings the gene or drug class. They could have differnt headers, so we work on column number
	#We set the first colum as row names and then remove it
	rownames(countdata)<-countdata[,1]
	countdata<-countdata[,-1]
	#Heavy approximation! DESeq2 wants integers, and we have to round to integers
	countdata<-round(countdata,0)
	metadata<-data.frame(colnames(countdata))
	rownames(metadata)<-metadata[,1]
	names(metadata)<-"Soil"
	#We only shorten the output names
 #browser()
	metadata$Soil<-substr(metadata$Soil,1,25)
 
	if(sum(rownames(metadata)==names(countdata))!=ncol(countdata)) stop("Names in countdata do not match names in metadata")
	ddsHTSeq <- DESeqDataSetFromMatrix(countData  = countdata,
                                        colData = metadata,
                                        design= ~ Soil)
	#Remove low abundance drug classes
    keep <- rowSums(counts(ddsHTSeq)) >= 5 
    ddsHTSeq<-ddsHTSeq[keep,]
	#Since we have no replicates (because we merged them), we NEED to use "blind=TRUE"
	#We could also rethink and see if we want to input the ppm and DO NOT perform normalization.
	#Right now we keep the vsd. We don't want to use this data for DE, so we are not scared
 	vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=TRUE)
	pdf(paste(graphdir,prefname,"_PCA.pdf",sep=""))
    myplot<-plotPCA(vsd,intgroup="Soil")
    print(myplot)
	dev.off()
  somma=rowSums(assay(vsd))
  mapvsd=assay(vsd)
  mapvsd=mapvsd[order(somma, decreasing = T),][1:min(nrow(mapvsd),50),]
pheatmap(mapvsd, cluster_rows=FALSE, fontsize=5, cellwidth=6, cellheight=4, filename=paste0(graphdir,prefname,"_vsd_heatmap.pdf"))

  keep <- rowSums(counts(ddsHTSeq)) >= 5 
    ddsHTSeq<-ddsHTSeq[keep,]
    browser()
	ddsHTSeq<-DESeq(ddsHTSeq)
 browser()
	#We are now comparing all vs all
  	to.contrast<-unique(metadata$Soil)
	compare<-combn(to.contrast,2)
    for(aaa in 1:ncol(compare))
    {
 		res <- results(ddsHTSeq,contrast=c("Soil",compare[2,aaa],compare[1,aaa]))
        res<-res[order(res$padj),]
        resfile<-paste(outdir,condition,"_",compare[2,aaa],"_vs_",compare[1,aaa],".txt",sep="")
		write.table(res,resfile,sep="\t",quote=F)
		mysig<-sum(res$padj<=0.05,na.rm=T)
		cat(compare[2,aaa],"vs",compare[1,aaa],",",mysig,"significant genes out of",nrow(res),"\n")
	}
}
runDEseq()
