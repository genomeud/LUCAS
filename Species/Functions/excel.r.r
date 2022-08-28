library(openxlsx)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel")
test<-read.xlsx("abundant_F")
average_soil_type <- read.xlsx("average_soil_type.xlsx")
average_soil_type <- average_soil_type[1:10,]