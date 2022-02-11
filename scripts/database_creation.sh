DATABASE=/mnt/vol1/projects/LUCAS/test_run/database/homolog_model/nucleotide_fasta_protein_homolog_model.fasta
DATABASENAME=nucleotide_fasta_protein_homolog_model
#kraken, we need to add a dummy "taxid" to every database entry 
awk -F "|" '{if($1 ~/^>/) {print $1"|"$3"|kraken:taxid|32630"}else{print $0}}' $DATABASE > ${DATABASENAME}_kraken.fasta
kraken2-build --download-taxonomy --db $DATABASENAME
kraken2-build --add-to-library ${DATABASE%.*}_kraken.fasta --db $DATABASENAME
kraken2-build --build --threads 24 --db $DATABASENAME
#bwa index
bwa index -p $DATABASENAME -a bwtsw $DATABASE
#dustmasker, the output will be used by samtools to mask low-abundance regions from BAM files.
dustmasker -in $DATABASE -outfmt acclist | awk -F ' ' '{print $1 "\t" $(NF-1) "\t" $NF}' | sed 's/>//g' > default_nucleotide_fasta_protein_homolog_model.fasta.txt
# .fai version of the database. Input of the faidx_to_gff.r function
samtools faidx $DATABASE
