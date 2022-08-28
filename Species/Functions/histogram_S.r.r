library(openxlsx)
library(data.table)
library(tibble)
library(ggplot2)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#inseriamo un file con estensione .txt
species<-fread("ratio_percentage_S.txt", data.table=F)
#creo un data frame con colonne e righe inverse rispetto a species
g<-data.frame(t(species))
#elimino la prima riga
d<-g[-1,]
#elimino le colonne che non ci interessano
d$X1<- d$X2<- d$X3<-NULL
#indico la colonna X4 come numero 
d$X4<-as.numeric(d$X4)
#modifico il nome della colonna
colnames(d)<- "ratio_percentage_classified_tot"
#creo l'istogramma con x=classified_tot_percentage y=numero di campioni con quella x
graph<-ggplot(data = d, aes(x = ratio_percentage_classified_tot)) + geom_histogram(bins = 30, fill = "seagreen", color="black")+
theme_minimal() + 
labs(x='Classificati / Totale', y='Conteggio', title='Istogramma delle specie classificate')
ggsave("Istogramma_classificati.pdf", width=30,units="cm",  plot=graph, device="pdf")