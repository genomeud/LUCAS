library(openxlsx)
library(tibble)
library(tidyr)
library(ggplot2)
library(strex)
#library(tidyverse)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel")
#caricato file di interesse
average_soil_type <- read.xlsx("average_soil_type.xlsx")
#salvataggio delle sole prime 10 righe
average_soil_type <- average_soil_type[1:10,]
k<-average_soil_type
k$Woodland<-k$Grassland<-k$Others<-k$Cropland<-NULL
colnames(k)<-"Name"
for ( i in 1 : nrow(k))
{
k[i,]<-str_before_nth(k[i,], "_", 2)
}
tibble_average<- average_soil_type %>%
pivot_longer(names_to = "Soil_type",
values_to = "percentage",
cols = -Name_Taxon_ID)
graph <- ggplot(tibble_average, aes(fill=Name_Taxon_ID, y=percentage, x=Soil_type)) + 
  geom_bar(position='stack', stat='identity') +
   theme_minimal() + 
  labs(x='Tipologia di suolo', y='Percentuale', fill='Nome', title='Specie più abbondanti per tipo di suolo') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_discrete( labels = k$Name)
  
ggsave("Specie_più_abbondanti_per_tipo_di_suolo.pdf", width=30,units="cm",  plot=graph, device="pdf")