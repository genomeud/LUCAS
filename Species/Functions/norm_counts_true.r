#library(pheatmap)
library(DESeq2)
library(data.table)
library(openxlsx)
library(gtools)
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
species <- fread("reduced_raw_counts_S.txt",data.table=F)
rownames(species)<- paste(species$Name, species$Taxon_ID, sep= "_")
species$Name<- species$Taxon_ID<- species$abundance<-NULL
species[is.na(species)]=0
species<- species[2:nrow(species),]

species <- species[, mixedsort(names(species))]
metadata<-read.xlsx("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel/metadata_plants.xlsx")
metadata$BARCODE_ID <- paste0("L", metadata$BARCODE_ID)
metadata<- metadata[,c ("BARCODE_ID", "Code")]
names(metadata)[2]<-"condition"
species<- species [, names(species) %in% metadata$BARCODE_ID]
rownames(metadata)<- metadata$BARCODE_ID
if(sum(rownames(metadata) ==names(species))!= ncol(species)) print("NO")
ddsHTSeqq <- DESeqDataSetFromMatrix(countData  = species,
                                        colData = metadata,
                                        design = ~ condition)
ddsHTSeqq <- DESeq(ddsHTSeqq)							
write.table(counts(ddsHTSeqq, normalized=TRUE), "norm_counts_true.txt", quote=F, sep="\t")

