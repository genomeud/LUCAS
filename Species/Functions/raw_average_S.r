library(data.table)
library(openxlsx)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
species <- fread("reduced_raw_counts_S.txt",data.table=F)
#elimono la riga unclassified
species<- species[2:nrow(species),]
#unisco specie e Taxon_ID
rownames(species)<- paste(species$Name, species$Taxon_ID, sep= "_")
#cancello le due colonne appena unire in rownames
species$Name <- species$Taxon_ID <- NULL
#csostituisco tutti gli NA con 0 per evitare errori 
species[is.na(species)]=0
metadata<-read.xlsx("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel/metadata_plants.xlsx")
metadata$BARCODE_ID <-paste0 ("L", metadata$BARCODE_ID)
#selezioniamo e sovrascriviamo solo le due colonne di interesse BARCODE_ID e LC_simpl_2018
metadata<- metadata[, c("BARCODE_ID", "LC_simpl_2018")]
#modifichiamo il nome della seconda colonna
names(metadata)[2]<- "condition"
#salvataggio di controllo
final <- species
#associo al vettore contition solo i valori singoli presenti all'interno della colonna condition
list_condition<- unique(metadata$condition)
#con il ciclo inserisco i
for ( i in 1 : length(list_condition))
{
subsetted<- subset(metadata, metadata$condition == list_condition [i])
middle <- species[, names(species) %in% subsetted$BARCODE_ID]
final[,list_condition[i]]<- apply ( middle, 1, mean)
}
#salva solo i dati di interesse
final<-final[,names(final) %in% list_condition]
final<-cbind(Name_Taxon_ID=rownames(final),final)
write.table(final, "raw_average_S.txt", row.names=F, quote=F)