#install.package("tibble")
library(tibble)
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
install.packages('bit64')
f630<-fread("reduced_raw_counts_F.txt",data.table = F)
summary(f630)
f630[is.na(f630)]=0
#sostituisco tutti i NA con 0 per avitare errori con l'utilizzo della funzione summary
F_summ<-as.data.frame(apply(f630[3:length(colnames(f630))],2,summary))
#converto il risultato della funzione summary effettuato dalla terza colonna all'ultima di f630 in un data frame con la funzione as.data.frame, per poter esportare su exel il risultato.
row.names(F_summ)<-gsub(" ","_",row.names(F_summ))
#nella prima colonnna di F_summ cambio tutti i caratteri " " in "_" 
F_summ<- tibble::rownames_to_column(F_summ, "Stat")
#attribuisco alla prima colonna il nome "Stat"
write.table(F_summ, "stat_summary_F.txt",row.names=F,quote=F)
test<-sapply(f630,function(x) head(f630$Name [order(x, decreasing=TRUE)],10))