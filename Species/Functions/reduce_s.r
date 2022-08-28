setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
library(data.table)
species <- fread("raw_counts_S.txt",data.table=F)
species$abundance<-rowSums(species[3:length(colnames(species))],na.rm=T)
s630<-subset(species,abundance>629)
s630<-s630[order(-s630$abundance),]
write.table(s630, "reduced_raw_counts_S.txt",row.names=F, quote=F)