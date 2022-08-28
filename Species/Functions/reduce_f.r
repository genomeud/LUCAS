library(data.table)
family<-fread("raw_counts_F.txt",data.table=F)
#attribuisco alla variabile family il contenuto del file raw_counts_F
family$abundance<-rowSums(family[3:length(colnames(family))],na.rm=T)
#creo una nuova colonna in family chiamada abundance, dove faccio la somma delle righe a partire dalla colonna tre all'ultima colonna di family. Per evitare problemi con l'uso della funzione rowSums indichiamo na.rm = T 
f630<-subset(family,abundance>629)
#tutte le righe in cui il valore contenuto nella colonna abundance>629 vengono salvate in f630
f630<-f630[order(-f630$abundance),]
#ordino f630 in funzione decrescente della funzione abundance
write.table(f630, "reduced_raw_counts_F.txt",row.names=F, quote=F)
#salvo il risultato in reduced_raw_counts_F. la funzione row.names=F permette di non salvare la colonna contenete i nomi delle righe. quote=F permette di salvare i caratteri non compresi fra le virgolette.