library(data.table)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file reduced_raw_counts_F.txt nel programma
species <- fread("reduced_raw_counts_S.txt",data.table=F)
#elimono la riga unclassified
species<- species[2:nrow(species),]
#elimino la colonna abundance
species$abundance <- NULL
#sostituisco tutti gli NA 
species[is.na(species)]=0
species$L471 <- species$L593 <- NULL
b_species<- species
#inseriamo due file uno con estensione .txt
tot<-fread("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle/countreads.txt", data.table=F)
#creo un vettore con il stesso numero di righe di species
ratio <- species$L1 
for (j in 3: ncol(species))
{
for (i in 1 : nrow(species))
{
ratio[i] <- (species[i,j] / tot$V2[j-2])*100
}
species[,j]<- ratio 
}
write.table(species, "percentage_counts.txt",row.names=F, quote=F)