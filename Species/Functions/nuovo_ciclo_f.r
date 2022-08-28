#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA'/Tesi/Analisi/Tabelle")
#caricamento della libreria per facilitare l'interazione con il formato del file
library(data.table)
#caricamento del file raw_counts_F.txt nel programma
family_reduced<-fread("raw_counts_F.txt",data.table=F)
#visualizzo le prime 4 righe e le prime 4 colonne
family_reduced[1:4,1:4]
#considerato che nella colonna Name ci sono dei doppioni, unisco nella stessa la colonna 1 e la colonna 2 (contentente il taxon_ID) all'interno della colonna Name
family_reduced$Name <- paste( family_reduced[,1] , family_reduced[,2], sep = "_") 
#elimono la colonna Taxon_ID
family_reduced$Taxon_ID <- NULL
#modifichiamo il nome delle righe rendendolo uguale alla colonna creata 1
row.names(family_reduced)<- family_reduced[,1]
#elimono la colonna 1 Name_ID
family_reduced$Name <- NULL
#sovrascrivo il data.frame per non prendere in considerazione il nome della righa contenente "unclassified_0"
family_reduced <- family_reduced [!(row.names(family_reduced) %in% c("unclassified_0")),]
#sostiuiamo tutti i NA dal data.frame family_reduced con 0
family_reduced[is.na(family_reduced)]=0
#creo un duplicato del data.frame chiamandolo final
final<-family_reduced
#utilizzando un ciclo ordiniamo una colonna alla volta, in ordine decrescente, e nel data.frame final sostituiamo la colonna corrispondente con il row.names
for (i in 1:ncol (family_reduced))
{
family_reduced <- family_reduced [order(family_reduced[,i], decreasing = TRUE) ,]
final[,i] <- rownames( paste( rownames (family_reduced) , family_reduced[,i], sep = "_"))
}
#salviamo senza considereare i row.names e togliendo le " con quote
write.table(final, "abundant_F.txt", row.names=F, quote=F)