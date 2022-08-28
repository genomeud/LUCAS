# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/analysis/all/reduced_0_stats.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-M", "--metafile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/metadata_plants.xlsx", 
              help="List of genes close to Vandal elements [default= %default]", metavar="character"),
  make_option(c("-R", "--readfile"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/metadata/countreads.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/AMR_complete/", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-C", "--condition"), type="character", default="Code,LC_simpl_2018,Country,Biogeographic_regions", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-G", "--graphdir"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/DESeq2/AMR_complete/", 
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
  library("stringr")
  library("dplyr")
	metadata<-read.xlsx(metafile)
  #metadata$LC1_2018 = substr(metadata$LC1_2018,1,nchar(metadata$LC1_2018)-1)
	countdata<-fread(abundance,data.table=F)
  rownames(countdata)<-countdata$Name
  countdata$Drug_class=countdata$Antibiotic=countdata$Resistance_mechanism=countdata$Name=NULL
	#rownames(countdata)<-countdata$Drug_class
	#countdata$Drug_class<-NULL
	metadata$BARCODE_ID<-paste0("L",metadata$BARCODE_ID)
	#Condition based on soil type: LC_simpl_2018
	#No covariates
	#We store bigmeta in memory for later re-use
	bigmeta<-metadata
allcond<-unlist(strsplit(condition,","))
options(scipen = 999)
	for(bbb in 1:length(allcond))
	{
	#if(bbb==2) browser()
	metadata<-bigmeta
	condition<-allcond[bbb]
	dir.create(paste0(outdir,condition))
	mydir<-paste0(outdir,condition,"/")
 padj="/padj"
  pdir<-paste0(outdir,condition,padj,"/")
  dir.create(pdir)
	metadata<-metadata[,c("BARCODE_ID",condition)]
	names(metadata)[2]<-"condition"
	readcount<-fread(readfile,data.table=F)
 #keep only the samples present in metadata, this allow to input, for instance, only the samples tagged as "italy" in the metadata and perform the analysis only on those samples
  readcount<-readcount[readcount$V1%in%metadata$BARCODE_ID,]
	countdata<-countdata[,names(countdata)%in%metadata$BARCODE_ID]
	metadata<-metadata[metadata$BARCODE_ID%in%metadata$BARCODE_ID,]
	rownames(metadata)<-metadata$BARCODE_ID
	#Add (unmapped) read counts to the count table (needed to normalize)
	#Trick to sort the read counts in the same order of the gene counts based on sample names
	mc<-match(names(countdata),readcount[,1])
	sreadcount<-readcount[mc,]
	tcount<-sreadcount[,2]
	totmapped<-apply(countdata,2,sum)
	unmapped<-tcount-totmapped
	unmapped_countdata<-data.frame(rbind(countdata,unmapped),stringsAsFactors=F)
	rownames(unmapped_countdata)[nrow(unmapped_countdata)]<-"Unmapped"
	#Check if metadata has NA and in case assign them to Unknown
	metadata$condition[is.na(metadata$condition)]<-"Unknown"
	#Strip leading and trailing whitespace from condition
	metadata$condition<-gsub("/","",metadata$condition)
	metadata$condition<-gsub(" ","",metadata$condition)
	#Check that samples are exactly the same and in the same order between counts and metadata
	if(sum(rownames(metadata)==names(countdata))!=ncol(countdata)) stop("Names in countdata do not match names in metadata")
	ddsHTSeq <- DESeqDataSetFromMatrix(countData  = countdata,
                                        colData = metadata,
                                        design= ~ condition)
	#Remove low abundance drug classes
    keep <- rowSums(counts(ddsHTSeq)) >= 5 
    ddsHTSeq<-ddsHTSeq[keep,]
    #cat("starting estimateSizeFactors \n")
    #pino<-estimateSizeFactors(ddsHTSeq, type = "iterate")
    ddsHTSeq<-DESeq(ddsHTSeq)
    write.table(counts(ddsHTSeq,normalized=TRUE),paste(mydir,"norm_counts.txt",sep=""),sep="\t",quote=F,col.names=NA)
    write.table(counts(ddsHTSeq,normalized=FALSE),paste(mydir,"raw_counts.txt",sep=""),sep="\t",quote=F,col.names=NA)
 unmapped_ddsHTSeq <- DESeqDataSetFromMatrix(countData  = unmapped_countdata,
                                        colData = metadata,
                                        design= ~ condition)
    ppm<-counts(unmapped_ddsHTSeq,normalized=FALSE)
	ppm<-ppm*1000000
	ppm<-t(apply(ppm,1,"/",tcount))
	write.table(ppm,paste(mydir,"ppm.txt",sep=""),sep="\t",quote=F,col.names=NA)
	vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
    # write.table(assay(vsd),paste(outdir,"VST.txt",sep=""),sep="\t",quote=F,col.names=NA)
	png(paste(mydir,"PCA.png",sep=""),height=10,width=10,units="cm",res=600,type="cairo")
    myplot<-plotPCA(vsd,intgroup="condition")
    print(myplot)
	dev.off()
	
	#Perform my own PCA analysis to retrieve the proportion of variance and write it in the axis
	mydat<-data.frame(assay(vsd),check.names=F)
	myvar<-apply(mydat,1,var)
	mydat<-mydat[order(myvar,decreasing=T),]
	#mydat<-mydat[1:500,]
	mypc<-prcomp(t(mydat))
	percentVar <- round(100*mypc$sdev^2/sum(mypc$sdev^2),0)
	myxlab<-paste("PC1: ",percentVar[1],"% variance",sep="")
	myylab<-paste("PC2: ",percentVar[2],"% variance",sep="")

	pdf(paste0(mydir,"VST_PCA.pdf"))
	#plotPCA(vsd, intgroup=c("condition"))
	plca<-plotPCA(vsd, intgroup=c("condition"),returnData=T)
 	graph<-ggplot(plca, aes(x=PC1, y=PC2, color=group),size=3)+geom_point(size=3) +geom_text_repel(aes(label=name),size=4) + coord_fixed() + theme(legend.text=element_text(size=10)) +xlab(myxlab) + ylab(myylab) 
	print(graph)
	dev.off()
	# pdf("dispersion_fit.pdf")
	# plotDispEsts(dd_1)
	# dev.off()

	#Build heatmap 
	newcounts<-counts(ddsHTSeq,normalized=TRUE)
	#We don't want to use unmapped, and we set them to zero
	newcounts=newcounts[rownames(newcounts)!="Unmapped",]
	keepme <- order(rowMeans(newcounts), decreasing=TRUE)[1:min(nrow(newcounts),50)]
	pcond<-data.frame(condition=metadata$condition,row.names=rownames(metadata))
	#pheatmap(newcounts[keepme,], cluster_rows=FALSE, annotation_col = pcond, fontsize=4, cellwidth=6, cellheight=4, filename=paste0(graphdir,"50-most-abundant-genes_clust.pdf"))
  #Build heatmap on vsd-corrected data
  somma=rowSums(assay(vsd))
  mapvsd=assay(vsd)
  mapvsd=mapvsd[order(somma, decreasing = T),][1:min(nrow(mapvsd),50),]
    pcond_color=pcond
  pcond_color$condition=as.factor(pcond$condition)
  ann_colors <- list(condition = c(Cropland = "Gold3",Grassland = "Green",Woodland = "Brown"))
  
pheatmap(mapvsd, cluster_rows=FALSE, annotation_col = pcond, fontsize=4, cellwidth=6, cellheight=4, filename=paste0(mydir,"vsd_50-most-abundant-genes_clust.pdf"))

    #dev.off()
	#Plot distance matrix as heatmap
	sampleDist <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDist)
	#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
	rownames(sampleDistMatrix) <- colnames(vsd)
	colnames(sampleDistMatrix) <- colnames(vsd)
	colors <- colorRampPalette ( rev( brewer.pal(9 , "Blues")) )( 255)
	pheatmap(sampleDistMatrix, clustering_distance_rows =sampleDist, annotation_col = pcond, annotation_row = pcond, clustering_distance_cols =sampleDist, filename=paste0(mydir,"VST_samples-distances.pdf"))
	#dev.off()

	#We are now comparing all vs all
	to.contrast<-unique(metadata$condition)
	compare<-combn(to.contrast,2)
    for(aaa in 1:ncol(compare))
    {
 		res <- results(ddsHTSeq,contrast=c("condition",compare[2,aaa],compare[1,aaa]))
        res<-res[order(res$padj),]
        pres<-subset(res, padj < 0.05)
        res=data.frame(class=row.names(res),res)
        pres=data.frame(class=row.names(pres),pres)
        resfile<-paste(mydir,compare[2,aaa],"_vs_",compare[1,aaa],".txt",sep="")
        presfile<-paste(pdir,compare[2,aaa],"_vs_",compare[1,aaa],".txt",sep="")
		write.table(res,resfile,sep="\t",quote=F,row.names=F)
   write.table(pres,presfile,sep="\t",quote=F,row.names=F)
		mysig<-sum(res$padj<=0.05,na.rm=T)
		cat(compare[2,aaa],"vs",compare[1,aaa],",",mysig,"significant genes out of",nrow(res),"\n")
   
	}
}
}
runDEseq()
