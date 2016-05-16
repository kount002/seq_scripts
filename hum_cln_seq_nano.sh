#processing of MiSeq/Nano raw data for human samples
#modified from hum_cln_seq.sh

#part1. preprocessing QC: fastqc

PDIR=$(pwd)

mkdir -p /nfs/gems_sata/tedder/evgueni/anal/human_processing
[ -e Processed_data ] || ln -s /nfs/gems_sata/tedder/evgueni/anal/human_processing/ Processed_data

mkdir -p "$PDIR"/qc/qc_pre

(/home/josh/collabs/software/FastQC/fastqc -casava -t 10 -o qc/qc_pre --noextract \
rawdata/*/*.fastq &> qc/qc_pre/fastqc.log)&
wait
echo 'pre-qc fastqc is complete'

#part2. adapter removal: cutadapt alone or wt combination with Trim Galore

#run cutadapt

for i in rawdata/* ; do
si=$(echo $i | cut -d "/" -f2)

mkdir -p "$PDIR"/Processed_data/"$si"

(
file1="$i"/"$si"_R1.fastq
file2="$i"/"$si"_R2.fastq 

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

rm "$PDIR"/Processed_data/"$si"/"$si"_atmp1.fastq.gz "$PDIR"/Processed_data/"$si"/"$si"_atmp2.fastq.gz

echo 'Adapters are removed'

#part3. post proceessing QC: fastqc

mkdir -p qc/qc_post

(/home/josh/collabs/software/FastQC/fastqc -t 10 -o qc/qc_post --noextract \
"$PDIR"/Processed_data/???N/*.fastq.gz &> qc/qc_post/postqc.log)&

echo 'Post-qc id has started'

#part4. tophat assembly and htscounts

for i in Processed_data/???N
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
echo tophat is done

for i in Processed_data/???N
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
echo htscounts are complete!

#part5. joining htscount output files into single file in preparation for R processing

for i in Processed_data/???N ; do
si=$(echo $i | cut -d "/" -f2)
awk '{print $1}' "$i"/"$si"_htscounts.txt > hts_count_merged_nano.txt
done

for i in Processed_data/???N ; do
si=$(echo $i | cut -d "/" -f2)
join hts_count_merged_nano.txt "$i"/"$si"_htscounts.txt >tmp_merged.txt
mv tmp_merged.txt hts_count_merged_nano.txt

done


#part6. DE Seq analysis
#part7. predominant and shared functional groups
#part8. prediminant and shared signaling networks
