#!/bin/bash
# give authorization
chmod 755 16S_workflow.sh

# assign input and output files
reads=$1
output=$2

# create the output folder
mkdir $2

# create the reads files list
find $1 -name "*.fastq.gz"| gunzip 
fastq_files=$(find $1 -name "*.fastq")

## create folder for fasta quality scores
mkdir $2/fastqc
for FILE in $fastq_files
do
	fastqc $FILE --outdir $2/fastqc
done

# create a dir for trimmed reads
mkdir $2/trimming

# perform trimming on reads
for elt in $(ls $1*.fastq |sed -r 's/_R[12].fastq//'|uniq)
do 
	file=$(basename $elt)
	java -jar ../tp_1/soft/AlienTrimmer.jar -if ${elt}_R1.fastq -ir ${elt}_R2.fastq -c ../tp_1/databases/contaminants.fasta -of $2/trimming/${file}_R1.trim.fastq -or $2/trimming/${file}_R2.trim.fastq -os $2/trimming/${file}_single.trim.fastq
done

# the next step: Vsearch
mkdir $2/vsearch
for elt in $(ls $2/trimming/*.trim.fastq |sed -r 's/_R[12].trim.fastq//'|uniq)
do

file=$(basename $elt)
suffixe=";sample=${file}" 
	../tp_1/soft/vsearch --fastq_mergepairs ${elt}_R1.trim.fastq --reverse ${elt}_R2.trim.fastq --label_suffix $suffixe --fastaout $2/vsearch/$file.merged.fasta
done

cat $2/vsearch/*.merged.fasta > $2/vsearch/amplicon.fasta

sed "s/ /_/g" $2/vsearch/amplicon.fasta > $2/vsearch/amplicon.clean.fasta

../tp_1/soft/vsearch --derep_fulllength $2/vsearch/amplicon.clean.fasta -sizeout --minuniquesize 10 --output $2/vsearch/amplicon.noSingleton.fasta
../tp_1/soft/vsearch --uchime_denovo $2/vsearch/amplicon.noSingleton.fasta --nonchimeras $2/vsearch/results.fasta
OTU_n="OTU_"
../tp_1/soft/vsearch --cluster_size $2/vsearch/results.fasta --id 0.97 --centroids $2/vsearch/centroids.fasta --uc $2/vsearch/clusters.uc --relabel $OTU_n

../tp_1/soft/vsearch --usearch_global $2/vsearch/results.fasta --otutabout $2/vsearch/amplicon_abun.txt --db $2/vsearch/centroids.fasta  --id 0.97

../tp_1/soft/vsearch --usearch_global $2/vsearch/centroids.fasta --userout  $2/vsearch/amplicon_annot.txt --db ../tp_1/databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields query+target

sed '1iOTU\tAnnotation' -i $2/vsearch/amplicon_annot.txt

