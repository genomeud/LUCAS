#Regroup the ARO in the output of trimming_to_htseq script by the drug_class targeted
Rscript ARO_to_class.r
#Produce the raw,normalized and ppm counts, plus Deseq2 analysis, PCA and heatmaps of "italy" tagged samples
Rscript DESeq2.r
#Group the samples in "italy" based on terrain type
Rscript terrain_type.r --token="cropland"
Rscript terrain_type.r --token="woodland"
Rscript terrain_type.r --token="grassland"
#Calculate the average ppm among samples of every drug class and join the terrain on a single table
for file in $(ls *land); do awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF-1); print$0 " "  sum}' $file > outfile; awk '{print$1 " " $NF}' outfile |  awk '{if (NR!=1) {print}}' > ${file}_average.txt; done
join cropland_average.txt woodland_average.txt > tmp
join tmp grassland_average.txt > tmp2
#Calculate the average of every row and select the 10 most abundant drug classes
awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF-1); print$0 " "  sum}' tmp2 | sort -n -k 4 -r | awk '{print $1 " " $2 " " $3 " " $4} NR==10{exit}' > ppm_average_terrain_type.txt
sed -i '1 i\Drug_Class Cropland Woodland Grassland' ppm_average_terrain_type.txt
#ppm_average_terrain_type.txt was exported and converted to ppm_average_terrain_type.xlsx 
#produce stacked bar plot of drug class in italy grouped by terrain type and plot the Log2FoldChange
Rscript bar_plot.r
