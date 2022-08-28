library(openxlsx)
library(data.table)
library(tibble)
library(ggplot2)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#inseriamo un file con estensione .txt
species<-fread("ratio_percentage_S.txt", data.table=F)
species$L471 <-species$L593<- NULL
#creo un data frame con colonne e righe inverse rispetto a species
g<-data.frame(t(species))
#elimino la prima riga
d<-g[-1,]
#elimino le colonne che non ci interessano
d$X1<- d$X2<- d$X4<-NULL
#indico la colonna X4 come numero 
d$X3<-as.numeric(d$X3)
#modifico il nome della colonna
colnames(d)<- "tot_reads"
#creo l'istogramma con x=classified_tot_percentage y=numero di campioni con quella x
graph<-ggplot(data = d, aes(x = tot_reads)) + geom_histogram(bins = 30, fill = "seagreen", color="black")+
theme_minimal() + 
labs(x='Numero di reads', y='Conteggio', title='Istogramma delle reads ottenute per campione')
ggsave("Istogramma_numero_di_reads.pdf", width=30,units="cm",  plot=graph, device="pdf")