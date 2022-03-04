#Regroup the ARO in the output of trimming_to_htseq script by the drug_class targeted
Rscript ARO_to_class.r
#Produce the raw,normalized and ppm counts, plus Deseq2 analysis, PCA and heatmaps of "italy" tagged samples
Rscript DESeq2.r
#Group the samples based on subclass
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Broadleaved_woodland"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Cereals"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Coniferous_woodland"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Dry_pulses_vegetables_and_flowers"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Fodder_crops"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Grassland_with_sparse_tree_shrub_cover"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Grassland_without_tree_shrub_cover"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Mixed_woodland"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Nonpermanent_industrial_crops"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Other_permanent_crops"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Permanent_crops_fruit_trees"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Root_crops"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Spontaneously_revegetated_surfaces"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Others"
#calculate the average
for file in $(ls); do awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF-1); print$0 " "  sum}' $file > outfile; awk '{print$1 " " $NF}' outfile |  awk '{if (NR!=1) {print}}' > average/${file%.*}_average.txt; done
cd average
#join the tables
join Cereals_average.txt Root_crops_average.txt > pingas
join pingas Nonpermanent_industrial_crops_average.txt > pino
join pino Dry_pulses_vegetables_and_flowers_average.txt > pingas
join pingas Fodder_crops_average.txt > pino
join pino Permanent_crops_fruit_trees_average.txt > pingas
join pingas Other_permanent_crops_average.txt > pino
join pino Broadleaved_woodland_average.txt > pingas
join pingas Coniferous_woodland_average.txt > pino
join pino Mixed_woodland_average.txt > pingas
join pingas Grassland_with_sparse_tree_shrub_cover_average.txt > pino
join pino Grassland_without_tree_shrub_cover_average.txt > pingas
join pingas Spontaneously_revegetated_surfaces_average.txt > pino
join pino Others_average.txt > bingas
#adding header
sed '1 i\ARO Name Drug_Class Antibiotic Resistance_Mechanism Cereals Root_crops Nonpermanent_industrial_crops Dry_pulses_vegetables_and_flowers Fodder_crops Permanent_crops_fruit_trees Other_permanent_crops Broadleaved_woodland Coniferous_woodland Mixed_woodland Grassland_with_sparse_tree_shrub_cover Grassland_without_tree_shrub_cover Spontaneously_revegetated_surfaces Others' bingas > average_raw_name_plants.txt
#ppm_average_terrain_type.txt was exported and converted to ppm_average_terrain_type.xlsx 
#same with soil type
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Woodland"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Cropland"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Others"
Rscript /mnt/vol1/projects/LUCAS/functions/select_test.r --token="Grassland"
for file in $(ls); do awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF-1); print$0 " "  sum}' $file > outfile; awk '{print$1 " " $NF}' outfile |  awk '{if (NR!=1) {print}}' > average/${file%.*}_average.txt; done
cd average
join Cropland_average.txt Grassland_average.txt > zino
join zino Others_average.txt > pingas
join pingas Woodland_average.txt > bingas
sed '1 i\ARO Name Drug_Class Antibiotic Resistance_Mechanism Cropland Grassland Others Woodland' bingas > average_drug_raw_soil.txt

#produce stacked bar plot of drug class grouped by terrain type and plot the Log2FoldChange
Rscript bar_plot.r
