library(pheatmap)
library(DESeq2)
library(data.table)
library(strex)
library(ggplot2)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
species <- fread("raw_average_S.txt",data.table=F)
#imposto il nome delle righe come il contenuto della colonna Name_Taxon_ID
rownames(species)<- species$Name_Taxon_ID
#elimino la ripetizione della colonna
species$Name_Taxon_ID<-NULL
#la funzione richiede solo numeri interi, per cui arrotondo i valori contenuti nel data.frame
species<-round(species, 0)
#creo un dataframe che contenga i tipi di suoli e ne modifico il nome delle colonne e delle righe
metadata<-data.frame(colnames(species))
rownames(metadata)<- metadata$colnames.species
names(metadata)<- "soil"
#trasformazione e normalizzazione dei dati necessria per poter visualizzare correttemente i grafici
ddsHTSeq <- DESeqDataSetFromMatrix(countData  = species,
                                        colData = metadata,
                                        design= ~ soil)
#normalizzazine estetica per la visualizzazione
vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=TRUE)				
mapvsd=assay(vsd)
somma<- rowSums(mapvsd)
mapvs<- mapvsd[order(somma, decreasing=T),][1:50,]


k<-rownames(mapvs)
for ( i in 1 : length(k))
{
k[i]<-str_before_last(k[i], "_")
}
rownames(mapvs)<- k

graph<-pheatmap(mapvs, main = "Cinquanta specie più abbondanti",fontsize=5, cellwidth=10, cellheight=4)

ggsave("Cinquanta_specie_più_abbondanti.pdf", width=20, units="cm",  plot=graph, device="pdf")