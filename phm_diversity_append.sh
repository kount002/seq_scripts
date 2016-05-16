# data files obtained form three crappy Micro MiSeq runs and were concatenated using append.run.sh

# The script evaluates quality, perform clean up (remove primer sequences at the end), puts the clean reads via tophat/hts count pipeline to evaluate number of genes (1) and number of reads for each gene (2). 

# This scrip is based on  hiseq_trio/hum_cln_seq.sh, next phm_diversity.sh

# START SCRIPT from within the directory
############################# 

#part1. append existing raw data files from three experiments.

PDIR=$(pwd)

#make processing directory to place work files
mkdir -p /nfs/gems_sata/tedder/evgueni/anal/phm_micro_append
[ -e Processed_data ] || ln -s /nfs/gems_sata/tedder/evgueni/anal/phm_micro_append/ Processed_data


#append files from tree runs to one big file in the processing directory
mkdir -p "$PDIR"/Processed_data/tmp
#remove preexisting files
rm -f "$PDIR"/Processed_data/tmp/*.fastq

for i in /home/kount002/Seq_data/Human/Micro/* ; do

# extract the name of the run directory
si=$(basename $i)


#append cycles for R1 and R2
(for j in /home/kount002/Seq_data/Human/Micro/"$si"/*_R1_* ; do
sn=$(basename $j | cut -c 1-5)
zcat $j >> "$PDIR"/Processed_data/tmp/"$sn"_R1_app.fastq ; done)&


(for j in /home/kount002/Seq_data/Human/Micro/"$si"/*_R2_* ; do
sn=$(basename $j | cut -c 1-5)
zcat $j >> "$PDIR"/Processed_data/tmp/"$sn"_R2_app.fastq ; done)&
wait

done
wait
echo 'Raw sequence files are concatenated'

#compress appended files

for i in "$PDIR"/Processed_data/tmp/* ; do
(gzip $i)& 
done
wait
echo 'Raw sequence files are gziped'

#part2. adapter removal: cutadapt alone or wt combination with Trim Galore
#run cutadapt



PDIR=$(pwd)
for i in "$PDIR"/Processed_data/tmp/*_R1_app.fastq.gz ; do
si=$(basename $i | cut -c 1-5 ) 
#cuts the file name

#dir=$(echo $si | sed 's/_R._tmp.fastq.gz/_app/') 
##cuts the sample name

mkdir -p "$PDIR"/Processed_data/"$si"

(
file1=$(ls "$PDIR"/Processed_data/tmp/"$si"*R1*)
file2=$(ls "$PDIR"/Processed_data/tmp/"$si"*R2*)

# R1 reads use universal primers for sequencing
#for R1 look for complement of 3-end of insert + Index adapter: GTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ( insert sequence starts with CCATGGCCGCCGAGAAC edd to universal adapter in reverse when cleaning for R2)
#Need to remove second adapter at 5' end that is insert start seq
~/.local/bin/cutadapt -g CCATGGCCGCCGAGAAC -a TGTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$si"/"$si"_atmp1.fastq.gz -p "$PDIR"/Processed_data/"$si"/"$si"_atmp2.fastq.gz \
$file1 $file2 \
&> Processed_data/"$si"/"$si"_cutadapt.log

#for R2 look for reverse complement of insert start and Universal primer (GTTCTCGGCGGCCATGG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT) 
~/.local/bin/cutadapt -a TGTTCTCGGCGGCCATGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$si"/"$si"_R2_trimmed.fastq.gz -p "$PDIR"/Processed_data/"$si"/"$si"_R1_trimmed.fastq.gz \
"$PDIR"/Processed_data/"$si"/"$si"_atmp2.fastq.gz "$PDIR"/Processed_data/"$si"/"$si"_atmp1.fastq.gz \
&>> Processed_data/"$si"/"$si"_cutadapt.log

rm "$PDIR"/Processed_data/"$si"/"$si"_atmp1.fastq.gz "$PDIR"/Processed_data/"$si"/"$si"_atmp2.fastq.gz
)&

echo "$si" 'adapters are removed'
done
wait
echo 'All adapters are removed'


#part3. post proceessing QC: fastqc
PDIR=$(pwd)


mkdir -p qc/qc_post

(/home/josh/collabs/software/FastQC/fastqc -t 10 -o qc/qc_post --noextract \
"$PDIR"/Processed_data/*/*.fastq.gz &> qc/qc_post/postqc.log)&

echo 'Post-qc id has started'

#part4. tophat assembly and htscounts

for i in Processed_data/*
do
si=$(echo $i | cut -d "/" -f2)
file1="$i"/"$si"_R1_trimmed.fastq.gz
file2="$i"/"$si"_R2_trimmed.fastq.gz

(
mkdir -p Processed_data/"$si"/tophat
tophat --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 -p 4 -r 150 -o Processed_data/"$si"/tophat \
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
rm -f Processed_data/tmp/*.txt

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
[ -e "$i"/"$si"_htscounts.txt ] && awk '{print $1}' "$i"/"$si"_htscounts.txt > hts_count_merged.txt
done

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
#head -2 "$i"/"$si"_htscounts.txt
[ -e "$i"/"$si"_htscounts.txt ] && join hts_count_merged.txt "$i"/"$si"_htscounts.txt >tmp_merged.txt

mv tmp_merged.txt hts_count_merged.txt
done



