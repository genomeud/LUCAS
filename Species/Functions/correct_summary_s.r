library(openxlsx)
library(data.table)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel")
#inseriamo due file uno con estensione .xlsx e uno con estensione .txt
correct_summary_S<-read.xlsx("tot_summary_S.xlsx")
correct_tot<-fread("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle/countreads.txt", data.table=F)
#elimino i campioni L471 e L593 
correct_summary_S$L471 <- correct_summary_S$L593 <- NULL
#modifichiamo il nome delle righe rendendolo uguale alla colonna creata 1
row.names(correct_summary_S)<- correct_summary_S[,1]
#elimono la colonna 1 Name_ID
correct_summary_S$X1 <- NULL
ifelse(colnames (correct_summary_S) == correct_tot$V1, "YES", "NO")
correct_summary_S[3,] <-correct_tot$V2
for (i in 1:length(correct_summary_S))
{
correct_summary_S[4,i]<- as.numeric(correct_summary_S[2,i]) / as.numeric(correct_summary_S[3,i])
}
write.table(correct_summary_S, "correct_summary_S.txt", quote=F)