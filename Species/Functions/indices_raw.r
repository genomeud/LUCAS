library(data.table)
library(vegan)
library(openxlsx)
library(tidyr)
library(ggplot2)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
norm_species<- fread("reduced_raw_counts_S.txt",data.table=F)
#cambio i rownames
rownames(norm_species)<- norm_species$Name
norm_species[is.na(norm_species)]=0
norm_species$L471<- NULL
norm_species$L593<- NULL
#cancello la colonna che contiene i dati salvati come rownames
norm_species$Name<- NULL
norm_species$Taxon_ID <- NULL
norm_species$abundance <- NULL
#creo un data frame con colonne e righe invertite rispetto a species
norm_species<-data.frame(t(norm_species))
norm_species$unclassified <- NULL
total<-rowSums (norm_species)
summ <-summary(total) [4]
soglia<- summ/1000
tot_species<- colSums (norm_species)
norm_species<- norm_species [, tot_species>soglia]
shannon<-diversity(norm_species, index = "shannon", MARGIN = 1, base = exp(1))
simpson<-diversity(norm_species, index = "simpson", MARGIN = 1, base = exp(1))
norm_species<-round(norm_species)
chao<-estimateR(norm_species)
write.table(shannon, "shannon.txt", quote=F)
write.table(simpson, "simpson.txt", quote=F)
write.table(chao, "chao.txt", quote=F)

metadata<-read.xlsx("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel/metadata_plants.xlsx")
metadata$BARCODE_ID <-paste0 ("L", metadata$BARCODE_ID)
shannon_frame<-data.frame (shannon)
shannon_frame <- tibble::rownames_to_column(shannon_frame, "BARCODE_ID")
shannon_meta <- merge (shannon_frame, metadata, by="BARCODE_ID")
shannon_wilcox <- pairwise.wilcox.test(fino$shannon, fino$LC_simpl_2018, p.adjust.method="none")
shannon_graph <- ggplot (data=shannon_meta, aes ( x = LC_simpl_2018, y = shannon)) + labs(y= "Shannon", x = "Tipologia di suolo") + ggtitle("Shannon box plot per tipologia di suolo") + geom_boxplot ()
ggsave("Shannon box plot.pdf", width=30,units="cm",  plot=shannon_graph, device="pdf")
shannon_graph <- ggplot (data=shannon_meta, aes ( , y = shannon)) + labs(y= "Shannon", x = "Conteggio") + ggtitle("Istogramma di Shannon") + geom_histogram () + coord_flip()
ggsave("Shannon histogram.pdf", width=30,units="cm",  plot=shannon_graph, device="pdf")

simpson_frame<-data.frame (simpson)
simpson_frame <- tibble::rownames_to_column(simpson_frame, "BARCODE_ID")
simpson_meta <- merge (simpson_frame, metadata, by="BARCODE_ID")
simpson_graph <- ggplot (data=simpson_meta, aes ( x = LC_simpl_2018, y = simpson))+labs(y= "Simpson", x = "Tipologia di suolo")+ ggtitle("Simpson box plot per tipologia di suolo") + geom_boxplot ()
ggsave("Simpson box plot.pdf", width=30,units="cm",  plot=simpson_graph, device="pdf")
simpson_graph <- ggplot (data=simpson_meta, aes ( , y = simpson)) + labs(y= "Simpson", x = "Conteggio") + ggtitle("Istogramma di Simpson") + geom_histogram () + coord_flip()
ggsave("Simpson histogram.pdf", width=30,units="cm",  plot=simpson_graph, device="pdf")

chao_frame<-data.frame (chao)
chao_frame<-chao_frame[ - c(1, 3, 4, 5),]
chao_frame<-data.frame(t(chao_frame))
chao_frame <- tibble::rownames_to_column(chao_frame, "BARCODE_ID")
chao_meta <- merge (chao_frame, metadata, by="BARCODE_ID")
chao_graph <- ggplot (data=chao_meta, aes ( x = LC_simpl_2018, y = S.chao1))+labs(y= "Chao", x = "Tipologia di suolo")+ ggtitle("Chao box plot per tipologia di suolo") + geom_boxplot ()
ggsave("Chao box plot.pdf", width=30,units="cm",  plot=chao_graph, device="pdf")
chao_graph <- ggplot (data=chao_meta, aes ( , y = S.chao1)) + labs(y= "Chao", x = "Conteggio") + ggtitle("Istogramma di Chao") + geom_histogram () + coord_flip()
ggsave("Chao histogram.pdf", width=30,units="cm",  plot=chao_graph, device="pdf")

chao_frame<-data.frame (chao)
chao_frame<-chao_frame[ - c(2, 3, 4, 5),]
chao_frame<-data.frame(t(chao_frame))
chao_frame <- tibble::rownames_to_column(chao_frame, "BARCODE_ID")
chao_meta <- merge (chao_frame, metadata, by="BARCODE_ID")
chao_graph <- ggplot (data=chao_meta, aes ( x = LC_simpl_2018, y = S.obs ))+labs(y= "Numero di specie", x = "Tipologia di suolo")+ ggtitle("Box plot ricchezza di specie per tipologia di suolo") + geom_boxplot ()
ggsave("Richness box plot.pdf", width=30,units="cm",  plot=chao_graph, device="pdf")
chao_graph <- ggplot (data=chao_meta, aes ( , y = S.obs )) + labs(y= "Ricchezza di specie", x = "Conteggio") + ggtitle("Istogramma della ricchezza di specie") + geom_histogram () + coord_flip()
ggsave("Richness histogram.pdf", width=30,units="cm",  plot=chao_graph, device="pdf")