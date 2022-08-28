#library(pheatmap)
library(DESeq2)
library(data.table)
library(openxlsx)
library(gtools)
#library(strex)
#library(ggplot2)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
species <- fread("reduced_raw_counts_S.txt",data.table=F)
#unisco specie e Taxon_ID e li uso come nomi delle righe
rownames(species)<- paste(species$Name, species$Taxon_ID, sep= "_")
#elimino la ripetizione della colonna
species$Name<- species$Taxon_ID<- species$abundance<-NULL
#elimono la riga unclassified
species<- species[2:nrow(species),]
#inseriamo un file uno con estensione .xlsx 
metadata<-read.xlsx("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel/metadata_plants.xlsx")
#aggiungo L alla colonna BARCODE_ID
metadata$BARCODE_ID <-paste0 ("L", metadata$BARCODE_ID)
#sovrascriviamo solo queste due colonne utili per effettuare la normalizzazione dei dati
metadata<- metadata[,c ("BARCODE_ID", "Code")]
species<- species [, names(species) %in% metadata$BARCODE_ID]
#metto BARCODE_ID come nome delle righe
rownames(metadata)<- metadata$BARCODE_ID
#contrrolliamo che ci sia lo stesso numero di colonne
if(sum(rownames(metadata) ==names(species))!= ncol(species)) print("NO")
#dà errore quindi usiamo una specifica funzione per metteli in ordine
names(species)<-mixedsort(names(species))
#tolgo tutti gli NA 
species[is.na(species)]=0
ddsHTSeq <- DESeqDataSetFromMatrix(countData  = species,
                                        colData = metadata,
                                        design = ~ Code)
#attribuisco alla variabile il risultato della funzione DESeq
ddsHTSeq <- DESeq(ddsHTSeq)
#salvo							
write.table(counts(ddsHTSeq, normalized=TRUE), "norm_counts_S.txt", quote=F, sep="\t")
