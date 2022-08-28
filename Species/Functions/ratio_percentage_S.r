library(openxlsx)
library(data.table)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel")
#inseriamo due file uno con estensione .xlsx
species<-read.xlsx("correct_summary_S.xlsx")
#attribusco a k la quarta riga di species, ratio
k<-species[4:4,]
#effettuo un ciclo per moltiplicare il valore a 100
for (i in 2:length(k))
{
k[,i]<- as.numeric(k[,i])*100
}
#modifico il nome della riga
k[1]<-"ratio_percentage"
#sovrascrivo alla quarta riga k
species[4:4,]<-k
write.table(species, "C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle/ratio_percentage_S.txt",row.names=F, quote=F)