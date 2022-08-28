options(scipen=999)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento della libreria per facilitare l'interazione con il formato del file
library(data.table)
#caricamento del file reduced_raw_counts_S.txt nel programma
species<-fread("reduced_raw_counts_S.txt",data.table=F)
species[is.na(species)]=0
#creo un vettore della stessa lunghezza di species
raw_tot<- colnames(species)
classified<- colnames(species)
unclassified<- species[1,]
#creo un vettore della stessa lunghezza di species
ratio<- colnames(species)
#attraverso un ciclo attribuisco al vettore la somma delle colonne di species
for (i in 3:length(species))
{
raw_tot[i]<-sum (species[,i])
}
#attraverso un ciclo effettuo il rapporto tra i non classificati e il totale
for (i in 3:length(raw_tot))
{
classified[i]<- as.numeric(raw_tot[i]) - species[1,i] 
}
for (i in 3:length(raw_tot))
{
ratio[i]<- species[1,i] / as.numeric(classified[i])
}
unclassified[ nrow(unclassified)+1,] <- classified
unclassified[ nrow(unclassified)+1,] <- raw_tot
unclassified[ nrow(unclassified)+1,] <- ratio
unclassified$Name <- unclassified$Taxon_ID <- NULL 
rownames(unclassified)<- c("unclassified" , "classified", "tot", "ratio")
write.table(unclassified, "tot_summary_S_bad.txt", quote=F)