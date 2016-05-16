# analyse number of sequence reads that are from custom sequencing primers using sequence set 4 

# data files obtained form Micro MiSeq run

# The script evaluates quality, perform clean up (remove primer sequences at the end), puts the clean reads via tophat/hts count pipeline to evaluate number of genes (1) and number of reads for each gene (2). 

# This scrip is based on  hiseq_trio/hum_cln_seq.sh

# START SCRIPT from within the directory
############################# 
# #part1. preprocessing QC: fastqc

#part2. adapter removal: cutadapt alone or wt combination with Trim Galore
#run cutadapt
PDIR=$(pwd)


for i in /home/kount002/Seq_data/Human/Micro/Phgm_diversity_3custom/*R1* ; do
si=$(echo $i | cut -d "/" -f8 ) 
#cuts the file name

dir=$(echo $si | sed 's/_S._L001//; s/_R._001.fastq.gz//') 
#cuts the sample name

mkdir -p "$PDIR"/Processed_data/"$dir"

(
file1=$(ls ~/Seq_data/Human/Micro/Phgm_diversity_3custom/"$dir"*R1*)
file2=$(ls ~/Seq_data/Human/Micro/Phgm_diversity_3custom/"$dir"*R2*)

# R1 reads use universal primers for sequencing
#for R1 look for complement of 3-end of insert + Index adapter: GTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ( insert sequence starts with CCATGGCCGCCGAGAAC edd to universal adapter in reverse when cleaning for R2)
#Need to remove second adapter at 5' end that is insert start seq
~/.local/bin/cutadapt -a TGTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz -p "$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz \
$file1 $file2 \
&> Processed_data/"$dir"/"$dir"_cutadapt.log

#for R2 look for reverse complement of insert start and Universal primer (GTTCTCGGCGGCCATGG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT) 
~/.local/bin/cutadapt -a TGTTCTCGGCGGCCATGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$dir"/"$dir"_R2_trimmed.fastq.gz -p "$PDIR"/Processed_data/"$dir"/"$dir"_R1_trimmed.fastq.gz \
"$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz \
&>> Processed_data/"$dir"/"$dir"_cutadapt.log

rm "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz "$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz
)&

done
wait


echo 'Adapters are removed'


#part4. tophat assembly and htscounts

for i in Processed_data/*
do
si=$(echo $i | cut -d "/" -f2)
file1="$i"/"$si"_R1_trimmed.fastq.gz
file2="$i"/"$si"_R2_trimmed.fastq.gz

(
mkdir -p Processed_data/"$si"/tophat
tophat --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 -p 8 -r 150 -o Processed_data/"$si"/tophat \
-G /nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
/nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \
$file1 $file2 &> Processed_data/"$si"/tophat/"$si"_tophat_stout.log
)&
done
wait
echo 'tophat is done'

for i in Processed_data/*
do
si=$(echo $i | cut -d "/" -f2)

samtools sort -on Processed_data/"$si"/tophat/accepted_hits.bam - | samtools view -h - | \
htseq-count -s no -m intersection-nonempty - \
/nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
> Processed_data/"$si"/"$si"_hts.txt

sample=$(echo 'sample_'$si)
sed "1s/.*/gene $sample/" "$i"/"$si"_hts.txt > Processed_data/"$si"/"$si"_htscounts.txt
#rm "$si"_hts.txt
done
echo 'htscounts are complete!'

#part5. joining htscount output files into single file in preparation for R processing

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
awk '{print $1}' "$i"/"$si"_htscounts.txt > hts_count_merged.txt
done

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
#head -2 "$i"/"$si"_htscounts.txt
join hts_count_merged.txt "$i"/"$si"_htscounts.txt >tmp_merged.txt
mv tmp_merged.txt hts_count_merged_3.txt

done


