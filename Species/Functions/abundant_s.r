library(data.table)
species_reduced<-fread("reduced_raw_counts_S.txt",data.table=F)
#considerato che nella colonna Name ci sono dei doppioni, unisco nella stessa la colonna 1 e la colonna 2 (contentente il taxon_ID) all'interno della colonna Name
species_reduced$Name <- paste( species_reduced[,1] , species_reduced[,2], sep = "_") 
species_reduced$Taxon_ID <- NULL
row.names(species_reduced)<- species_reduced[,1]
species_reduced$Name <- NULL
species_reduced <- species_reduced [!(row.names(species_reduced) %in% c("unclassified_0")),]
species_reduced[is.na(species_reduced)]=0
final<- species_reduced
for (i in 1:ncol(species_reduced))
{
species_reduced <- species_reduced[order(species_reduced[,i],decreasing = TRUE),]
final[,i] <- paste(rownames(species_reduced), species_reduced[,i], sep = "|")
}
vector_names<- colnames(final)
for (i in 1:length (vector_names))
{
vector_names[i]<- paste( vector_names[i], "|", vector_names[i], "_raw_count", sep="")
}
colnames (final)<- vector_names
final<- head (final, 50)
write.table(final, "abundant_S.txt", row.names=F, quote=F)