setwd(C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle)
library(data.table)
genus <- fread("raw_counts_G.txt",data.table=F)
genus$abundance<-rowSums(genus[3:length(colnames(genus))],na.rm=T)
g630<-subset(genus,abundance>629)
g630<-g630[order(-g630$abundance),]
write.table(g630, "reduced_raw_counts_G.txt",row.names=F, quote=F)