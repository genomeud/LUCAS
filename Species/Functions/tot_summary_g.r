options(scipen=999)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento della libreria per facilitare l'interazione con il formato del file
library(data.table)
#caricamento del file reduced_raw_counts_G.txt nel programma
genus<-fread("reduced_raw_counts_G.txt",data.table=F)
genus[is.na(genus)]=0
#creo un vettore della stessa lunghezza di genus
raw_tot<- colnames(genus)
classified<- colnames(genus)
unclassified<- genus[1,]
#creo un vettore della stessa lunghezza di genus
ratio<- colnames(genus)
#attraverso un ciclo attribuisco al vettore la somma delle colonne di genus
for (i in 3:length(genus))
{
raw_tot[i]<-sum (genus[,i])
}
#attraverso un ciclo effettuo il rapporto tra i non classificati e il totale
for (i in 3:length(raw_tot))
{
classified[i]<- as.numeric(raw_tot[i]) - genus[1,i] 
}
for (i in 3:length(raw_tot))
{
ratio[i]<- genus[1,i] / as.numeric(classified[i])
}
unclassified[ nrow(unclassified)+1,] <- classified
unclassified[ nrow(unclassified)+1,] <- raw_tot
unclassified[ nrow(unclassified)+1,] <- ratio
unclassified$Name <- unclassified$Taxon_ID <- NULL 
rownames(unclassified)<- c("unclassified" , "classified", "tot", "ratio")
write.table(unclassified, "tot_summary_G.txt", quote=F)