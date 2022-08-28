#install.package("tibble")
library(tibble)
setwd("C:/Users/utente/Desktop/UNIVERSITA'/Tesi/Analisi/Tabelle")
#install.packages('bit64')
s630<-fread("reduced_raw_counts_S.txt",data.table=F)
summary(s630)
s630[is.na(s630)]=0
S_summ<-as.data.frame(apply(s630[3:length(colnames(s630))],2,summary))
row.names(S_summ)<-gsub(" ","_",row.names(S_summ))
S_summ<- tibble::rownames_to_column(S_summ, "Stat")
write.table(S_summ, "stat_summary_S.txt",row.names=F,quote=F)
test<-sapply(s630,function(x) head(s630$Name [order(x, decreasing=TRUE)],10))