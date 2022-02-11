#Group each ARO by the drug class targeted
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--aro"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/analysis/all/5_stats",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="/mnt/vol1/projects/LUCAS/JRC_result/analysis/all/aro_to_class_tmp.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$aro)) {
  stop("WARNING: No aro specified with '-A' flag.")
} else {  cat ("aro is ", opt$aro, "\n")
  aro <- opt$aro  
  }


  if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output dir is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

#Function to assign reads mapping to each ARO to the correspnding drug class. 
#I one ARO matches to "n" drug classes, we attribute 1/n reads to each class.
aro_to_class<-function() 
{
    library(data.table)
	countdata<-fread(aro,data.table=F)
	uniq_class<-unique(unlist(strsplit(countdata$Drug_Class,";")))
	
	countdata$nclasses<-unlist(lapply(strsplit(countdata$Drug_Class,";"),length))
	#Remember: NEVER add column names starting with L!!!!!
	todiv<-names(countdata)[substr(names(countdata),1,1)=="L"]
	newdata<-countdata[,todiv]/countdata$nclasses
	rownames(newdata)<-countdata$ARO
	
	#Create data frame with drug classes and sample names
	drug_df<-matrix(0,nrow=length(uniq_class),ncol=length(todiv),dimnames=list(uniq_class,todiv))
	#Loop over AROs
	for(aaa in 1:nrow(newdata))
	{
	sclass<-unlist(strsplit(countdata$Drug_Class[aaa],";"))
	#Loop over antibiotic class for each ARO
		for(bbb in 1:length(sclass))
		{
		drug_df[sclass[bbb],]<-drug_df[sclass[bbb],]+unlist(newdata[aaa,])
		}
	}
	drug_df<-round(drug_df,0)
	newdrug_df<-data.frame(drug_df)
	newdrug_df$Drug_class<-row.names(newdrug_df)
	newdrug_df<-newdrug_df[,c("Drug_class",todiv)]
	write.table(newdrug_df,outfile,sep="\t",row.names=F,quote=F)
}
aro_to_class()
