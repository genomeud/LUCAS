options(scipen=999)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento della libreria per facilitare l'interazione con il formato del file
library(data.table)
#caricamento del file reduced_raw_counts_F.txt nel programma
family<-fread("reduced_raw_counts_F.txt",data.table=F)
family[is.na(family)]=0
#creo un vettore della stessa lunghezza di family
raw_tot<- colnames(family)
classified<- colnames(family)
unclassified<- family[1,]
#creo un vettore della stessa lunghezza di family
ratio<- colnames(family)
#attraverso un ciclo attribuisco al vettore la somma delle colonne di family
for (i in 3:length(family))
{
raw_tot[i]<-sum (family[,i])
}
#attraverso un ciclo effettuo il rapporto tra i non classificati e il totale
for (i in 3:length(raw_tot))
{
classified[i]<- as.numeric(raw_tot[i]) - family[1,i] 
}
for (i in 3:length(raw_tot))
{
ratio[i]<- family[1,i] / as.numeric(classified[i])
}
unclassified[ nrow(unclassified)+1,] <- classified
unclassified[ nrow(unclassified)+1,] <- raw_tot
unclassified[ nrow(unclassified)+1,] <- ratio
unclassified$Name <- unclassified$Taxon_ID <- NULL 
rownames(unclassified)<- c("unclassified" , "classified", "tot", "ratio")
write.table(unclassified, "tot_summary_F.txt", quote=F)