module load compilers/gcc/5.4.0_novabreed
module load aligners/bwa/0.7.13
module load aligners/blast/2.6.0+
#Function to produce the GFF file required by htseq
Rscript faidx_to_gff.r
DUSTMASKER=/mnt/vol1/projects/LUCAS/test_run/database/homolog_model/default_nucleotide_fasta_protein_homolog_model.fasta.txt;  
SAMPLEDIR=/mnt/vol1/projects/LUCAS/JRC_data
RESULTDIR=/mnt/vol1/projects/LUCAS/JRC_result
TABLEDIR=/mnt/vol1/projects/LUCAS/JRC_result/table
SWISSPROT=/mnt/vol1/biodb/ncbi/blastdb/2019-10-23/swissprot
KRAKEN_DATABASE=/mnt/vol1/projects/LUCAS/test_run/database/kraken/nucleotide_fasta_protein_homolog_model
KRAKEN=/mnt/vol1/projects/SRAome/kraken2/kraken2
BWADATABASE=/mnt/vol1/projects/LUCAS/test_run/database/homolog_model/nucleotide_fasta_protein_homolog_model.fasta
GFF=/mnt/vol1/projects/LUCAS/test_run/database/homolog_model/my_simple_gtf.gff
cd $SAMPLEDIR
for dir in $(ls -d */); do 
#trimming to remove adapters 
conda activate BBMap;  
for f in $(ls $dir*.fq.gz | sed -e 's/_1.fq.gz//' -e 's/_2.fq.gz//' | sort -u);do
bbduk.sh in1=${f}_1.fq.gz in2=${f}_2.fq.gz out1=$RESULTDIR/$dirx_${f}_1_trim.fq.gz out2=$RESULTDIR/$dirx_${f}_2_trim.fq.gz ref=${BBMAPDIR}/resources/adapters.fa ktrim=r ktrim=l qtrim=rl trimq=10;      
done;
done   
conda deactivate
# Merging reads obtained in different runs for the same sample 
cd $RESULTDIR
for dir in $(ls);  do   
count=$( ls $dir | wc -l)
if [ $count -eq 4 ]
then
if [ -e $dir/*H2MKGDSXY_L2_1* ]; then
cat $dir/*H2KT3DSXY_L3_1_trim.fq.gz $dir/*H2MKGDSXY_L2_1_trim.fq.gz > $dir/${dir}_1_merged.fq.gz
cat $dir/*H2KT3DSXY_L3_2_trim.fq.gz $dir/*H2MKGDSXY_L2_2_trim.fq.gz > $dir/${dir}_2_merged.fq.gz 
else
cat $dir/*32DSXY_L3_1_trim.fq.gz $dir/*KGDSXY_L3_1_trim.fq.gz > $dir/${dir}_1_merged.fq.gz
cat $dir/*32DSXY_L3_2_trim.fq.gz $dir/*KGDSXY_L3_2_trim.fq.gz > $dir/${dir}_2_merged.fq.gz  
fi 
fi 
done

for dir in $(ls);  do   
count=$( ls $dir | wc -l)    
if [ $count -eq 8 ]    
then   
if [ -e $dir/*H2MKGDSXY_L2_1* ]   
 then 
 cat $dir/*H2KT3DSXY_L3_1_trim.fq.gz  $dir/*H2MKGDSXY_L2_1_trim.fq.gz >  $dir/${dir}_1_merged.fq.gz 
 cat $dir/*H2KT3DSXY_L3_2_trim.fq.gz  $dir/*H2MKGDSXY_L2_2_trim.fq.gz >  $dir/${dir}_2_merged.fq.gz  
 cat $dir/*H2TVMDSXY_L2_1_trim.fq.gz $dir/*H2V7NDSXY_L2_1_trim.fq.gz >  $dir/meta9_${dir}_1_merged..3ta9_${dir}_2_merged.fq.gz 
 else 
 cat $dir/*32DSXY_L3_1_trim.fq.gz $dir/*KGDSXY_L3_1_trim.fq.gz > $dir/${dir}_1_merged.fq.gz  
 cat  $dir/*32DSXY_L3_2_trim.fq.gz $dir/*KGDSXY_L3_2_trim.fq.gz > $dir/${dir}_2_merged.fq.gz 
 cat  $dir/*H2TVMDSXY_L2_1_trim.fq.gz $dir/*H2V7NDSXY_L2_1_trim.fq.gz > $dir/meta9_${dir}_1_merged.fq.gz  
 cat  $dir/*H2TVMDSXY_L2_2_trim.fq.gz $dir/*H2V7NDSXY_L2_2_trim.fq.gz > $dir/meta9_${dir}_2_merged.fq.gz  
 fi  
 fi  
done
#First alignment against CARD using kraken
for dir in $(ls -d */);  do 
for f in $(ls $dir*_merged.fq.gz | sed -e 's/_1_merged.fq.gz//' -e 's/_2_merged.fq.gz//' | sort -u);  do   
$KRAKEN -t 32 --db $KRAKEN_DATABASE --paired --classified-out $dirx${f}_merged#.fq ${f}_1_merged.fq.gz ${f}_2_merged.fq.gz > ${f}.txt 
done
for dir in $(ls -d */);  do
for f in $(ls $dir*.fq | sed -e 's/_1.fq//' -e 's/_2.fq//' | sort -u); do
 #Second alignment against CARD using BWA
 bwa mem -t 24 $BWADATABASE $dirx${f}_1.fq $dirx${f}_2.fq | samtools view -bS - > $dirx${f}.bam;  
 done;  
#Masking of low complexity sequences
  for f in $(ls $dir*.bam); do 
  conda activate samtools
  samtools view $f  -b -h -o REMOVEME -U ${f%.*}_default.bam -L $DUSTMASKER
  done
  conda deactivate 
  #Reads count with htseq
   conda activate htseq
  for f in $(ls $dir*_default.bam); do 
 module load tools/samtools/0.1.19; 
 samtools sort -@ 8 -n $f ${f/.bam/_sortname}; 
 NAMEBAM=${f/.bam/_sortname}.bam; 
 htseq-count -a 0 -s no $NAMEBAM $GFF > ${f%.*}_htseq.txt; 
 done;
  #Editing to obtain the final tables
 for f in $(ls $dir*_htseq.txt); do
 awk -F '|' '{ print $5 " " $6}' $f | awk '{print $1 " " $NF}' | sort |sed '1,5d'> ${f%.*}_ID.txt
 join /mnt/vol1/projects/LUCAS/test_run/database/homolog_model/ID_novariant.txt ${f%.*}_ID.txt | sort -n -k3 -r > ${f%.*}_aro.txt
done
#bam to fasta and blast as quality control
for f in $(ls $dir*_default.bam); do
samtools view $f | awk '{print $3 "\n" $10}' |  sed '1~2s/^/>/' > ${f%.*}.fasta
conda deactivate
done
#blast
for f in $(ls $dir*_default.fasta); do
 blastx -db $SWISSPROT -query $f -num_threads 24 -out ${f%.*}_blast.txt  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -num_alignments 2
done
done
#Grouping all the htseq output into a single table and adding ARO code
 for dir in $(ls -d L*); do
  if [  -e $dir/meta*ID* ]; then
 cp $dir/meta*ID* $TABLEDIR
 else
 cp $dir/*ID* $TABLEDIR
 fi
 done
 cd $TABLEDIR
 for f in $(ls *ID*); do
 join id.txt $f > yo
 mv yo id.txt
 done
