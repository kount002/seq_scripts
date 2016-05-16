#processing of HiSeq raw data for human samples

#part1. preprocessing QC: fastqc

PDIR=$(pwd)

mkdir -p /nfs/gems_sata/tedder/evgueni/anal/human_processing
[ -e Processed_data ] || ln -s /nfs/gems_sata/tedder/evgueni/anal/human_processing/ Processed_data

mkdir -p "$PDIR"/qc/qc_pre

(/home/josh/collabs/software/FastQC/fastqc -casava -t 10 -o qc/qc_pre --noextract \
~/Seq_data/Human/HiSeq/*/*.fastq.gz &> qc/qc_pre/fastqc.log)&
wait
echo 'pre-qc fastqc is complete'

#part2. adapter removal: cutadapt alone or wt combination with Trim Galore

#append files
mkdir -p "$PDIR"/Processed_data/tmp

for i in /home/kount002/Seq_data/Human/HiSeq/* ; do
si=$(echo $i | cut -d "/" -f7)

[ -e  "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq ] && rm  "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq
[ -e  "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq ] && rm  "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq

(for j in /home/kount002/Seq_data/Human/HiSeq/"$si"/*_R1_* ; do
zcat $j >> "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq ; done)&

(for j in /home/kount002/Seq_data/Human/HiSeq/"$si"/*_R2_* ; do
zcat $j >> "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq ; done)&

done
wait
echo 'Raw sequence files are concatenated'

#run cutadapt

for i in /home/kount002/Seq_data/Human/HiSeq/* ; do
si=$(echo $i | cut -d "/" -f7)

mkdir -p "$PDIR"/Processed_data/"$si"

(
file1="$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq
file2="$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq 

#for R1 
~/.local/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q20 -m20 -e 0.1 \
-o Processed_data/"$si"/"$si"_atmp1.fastq.gz -p Processed_data/"$si"/"$si"_atmp2.fastq.gz \
$file1 $file2 \
&> Processed_data/"$si"/"$si"_cutadapt.log


#for R2
~/.local/bin/cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q20 -m20 -e 0.1 \
-o Processed_data/"$si"/"$si"_R2_trimmed.fastq.gz -p Processed_data/"$si"/"$si"_R1_trimmed.fastq.gz \
Processed_data/"$si"/"$si"_atmp2.fastq.gz Processed_data/"$si"/"$si"_atmp1.fastq.gz \
&>> Processed_data/"$si"/"$si"_cutadapt.log
)&

done 
wait

rm "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq \
"$PDIR"/Processed_data/"$si"/"$si"_atmp1.fastq "PDIR"/Processed_data/"$si"/"$si"_atmp2.fastq

rmdir "$PDIR"/Processed_data/tmp/
echo 'Adapters are removed'

#part3. post proceessing QC: fastqc

mkdir -p qc/qc_post

(/home/josh/collabs/software/FastQC/fastqc -t 10 -o qc/qc_post --noextract \
"$PDIR"/Processed_data/*/*.fastq.gz &> qc/qc_post/postqc.log)&

echo 'Post-qc id complete'

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

samtools sort -on Processed_data/"$si"/tophat/accepted_hits.bam - | samtools view -h - | \
htseq-count -s no -m intersection-nonempty - \
/nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
> Processed_data/"$si"/"$si"_hts.txt

sample=$(echo 'sample_' $si)
sed "1s/.*/gene $sample/" "$si"_hts.txt > Processed_data/"$si"/"$si"_htscounts.txt
#rm "$si"_hts.txt
)&

done
wait
echo tophat and htscounts are complete!

#part5. joining htscount output files into single file in preparation for R processing

touch hts_count_merged.txt

for i in Processed_data/*
si=$(echo $i | cut -d "/" -f2)

join hts_count_merged.txt "$i"/"si"_htscounts.txt >tmp_merged.txt
mv tmp_merged.txt hts_count_merged.txt

done
rm tmp_merged.txt 


#part6. DE Seq analysis
#part7. predominant and shared functional groups
#part8. prediminant and shared signaling networks
