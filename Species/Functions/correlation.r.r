#library(pheatmap)
library(DESeq2)
library(data.table)
library(openxlsx)
library(gtools)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
#posizionamento nella cartella corretta
setwd("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Tabelle")
#caricamento del file 
norm_species <- fread("norm_counts_true.txt",data.table=F)
#inseriamo un file uno con estensione .xlsx 
metadata<-read.xlsx("C:/Users/utente/Desktop/UNIVERSITA/Tesi/Analisi/Exel/metadata_plants.xlsx")
#aggiungo L alla colonna BARCODE_ID
metadata$BARCODE_ID <-paste0 ("L", metadata$BARCODE_ID)
#cambio i rownames
rownames(norm_species)<- norm_species$V1
#cancello la colonna che contiene i dati salvati come rownames
norm_species$V1<- NULL

tot<-rowSums(norm_species)
norm_species<-norm_species[order(tot,decreasing=T),][1:25,]
#creo un data frame con colonne e righe invertite rispetto a species
norm_species<-data.frame(t(norm_species))
#creo un vettore che abbia i rownames di norm_species
BARCODE_ID<-rownames(norm_species)
#inserico una colonna per prima che contiene il vettore BARCODE_ID
norm_species<-cbind (BARCODE_ID,norm_species)
#mi servono due vettori che contengano i nomi delle specie e i caratteri che voglio studiare (specie) (variabili)
pippo<-colnames(norm_species[2:ncol(norm_species)])
topolino<-c("GPS_LAT","GPS_LONG","Coarse_fragments","Clay_content","Sand_content","Silt_content","pH_CaCl2",
			"pH_H2O","Electrical_conductivity","Organic_carbon","Carbonate_content","Phosphorus_content","Total_nitrogen_content",
			"Extractable_potassium_content","Bulk_Density_0_10_cm","Bulk_Density_10_20_cm")
#merge che unisce due file in base alla colonna BARCODE_ID
tcount <- merge (norm_species, metadata, by="BARCODE_ID", sort=F)
#utilizziamo la funzione expanded.grid genera tutte le combinazioni possibili fra pippo e topolino
cdf<- data.frame (expand.grid(pippo, topolino),stringsAsFactors=F)
#rinomino le colonne
names(cdf) <- c ("Specie", "Proprieta")
#impongo che siano classificati come caratteri
cdf$Specie <- as.character (cdf$Specie)
cdf$Proprieta <- as.character (cdf$Proprieta)

for(aaa in 1:nrow(cdf))
	{
		pp<-cor.test(tcount[,cdf$Specie[aaa]],tcount[,cdf$Proprieta[aaa]],method="spearman")
  	cdf$rho[aaa]<-round(pp$estimate,3)
		cdf$pvalue[aaa]<-signif(pp$p.value,3)
	}
	
cdf$fdr<-p.adjust(cdf$pvalue,method="fdr")

#write.table(cdf, "correlation.txt", quote=F, row.names=F, sep="\t")

cormat<-corp<-matrix(NA,nrow=length(pippo),ncol=length(topolino),dimnames=list(pippo,topolino))
	#I run this second loop instead of using just one because I want to use fdr, and this is computed AFTER finishing the first loop
	for(aaa in 1:nrow(cdf))
	{
		cormat[cdf$Specie[aaa],cdf$Proprieta[aaa]]<-cdf$rho[aaa]
		corp[cdf$Specie[aaa],cdf$Proprieta[aaa]]<-cdf$fdr[aaa]
	}
	
sig<- cormat<=0.05
tsig<-apply (sig, 1, "sum")
cormat<-cormat[tsig>0,]
corp<-corp [tsig>0,]

large<-abs(cormat)>0.2
	tlarge<-apply(large,1,"sum")
	cormat<-cormat[tlarge>0,]
	corp<-corp[tlarge>0,]

col2 = colorRampPalette(c('darkred', 'white', 'blue3'))
pdf("megacorr_class.pdf", width=10)

corrplot(cormat, p.mat = corp, 
         tl.col="black", tl.cex=0.6,
         diag = TRUE, type = 'full',
         col=col2(256),
         insig = 'label_sig', sig.level = c(0.05), 
         title="Correlazione fra le venticinque specie più abbondanti e le proprietà del suolo",
		 mar=c(0,0.5,1.2,0.5),oma=c(0,0,0),
         pch.cex = 0.8, pch.col = 'black')
	dev.off()