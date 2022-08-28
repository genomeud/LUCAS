library(data.table)
library(vegan)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
norm_species<- fread("norm_counts_true.txt",data.table=F)
#cambio i rownames
rownames(norm_species)<- norm_species$V1
#cancello la colonna che contiene i dati salvati come rownames
norm_species$V1<- NULL
#creo un data frame con colonne e righe invertite rispetto a species
norm_species<-data.frame(t(norm_species))
shannon<-diversity(norm_species, index = "shannon", MARGIN = 1, base = exp(1))
simpson<-diversity(norm_species, index = "simpson", MARGIN = 1, base = exp(1))
norm_species<-round(norm_species)
chao<-estimateR(norm_species)