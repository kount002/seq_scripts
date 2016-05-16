# modified from phm_divesity_genelist.sh to analyse Miseq data for Severn Phage pipe libraries
# Provide a location for the data files

# modified from phm_diversity_append.sh, to work on single set of sequences from custom primer

# The script evaluates quality, perform clean up (remove primer sequences at the end), puts the clean reads via tophat/hts count pipeline to evaluate number of genes (1) and number of reads for each gene (2). 

# This scrip is based on  hiseq_trio/hum_cln_seq.sh, next phm_diversity.sh

# START SCRIPT from within the directory
############################# 


#part1. prepare the environment for the processing.
echo "The script_" $(basename $0) "_was started at:" $(date)

PDIR=$(pwd)

#make processing directory to place work files
mkdir -p /nfs/gems_sata/tedder/evgueni/anal/$(basename $PDIR)
[ -e Processed_data ] || ln -s /nfs/gems_sata/tedder/evgueni/anal/$(basename $PDIR)/ Processed_data

#part2. adapter removal: cutadapt alone or wt combination with Trim Galore
#run cutadapt

for i in /home/kount002/Seq_data/Human/MiPhage/*R1* ; do
si=$(basename $i) #cuts the file name
dir=$(echo $si | cut -d "_" -f1) #cuts the sample name
echo "Adapter processing of ....." "$dir", "$si"
mkdir -p "$PDIR"/Processed_data/"$dir"

(
#create input file variables
file1=$(ls ~/Seq_data/Human/MiPhage/"$dir"*R1*)
file2=$(ls ~/Seq_data/Human/MiPhage/"$dir"*R2*)

# R1 reads use universal primers for sequencing works with compressed data
#for R1 look for complement of 3-end of insert + Index adapter: GTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ( insert sequence starts with CCATGGCCGCCGAGAAC edd to universal adapter in reverse when cleaning for R2)
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

echo "$si" 'adapters are being removed'
done
wait

#check if the processed files are exitst, otherwise exit 2
ls Processed_data/*/*.fastq.gz &> /dev/null || ( echo "Files with removed adapter not found" && exit 2 )
echo 'All adapters were removed'


#part3. post proceessing QC: fastqc
mkdir -p qc/qc_post
(/home/josh/collabs/software/FastQC/fastqc -t 10 -o qc/qc_post --noextract \
"$PDIR"/Processed_data/*/*.fastq.gz &> qc/qc_post/postqc.log)&
echo 'Post-qc has started'

#part4. tophat assembly and htscounts

for i in Processed_data/*
do
si=$(basename $i) #cuts the sample/directory name

file1="$i"/"$si"_R1_trimmed.fastq.gz
file2="$i"/"$si"_R2_trimmed.fastq.gz
echo "Tophat processing ...." "$i", "$si"

(
mkdir -p Processed_data/"$si"/tophat
tophat --read-mismatches 4 --read-gap-length 2 --read-edit-dist 5 -p 4 -r 150 -o Processed_data/"$si"/tophat \
-G /nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
/nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \
$file1 $file2 &> Processed_data/"$si"/tophat/"$si"_tophat_stout.log
)&

echo "$si" "...is being processed by tophat" 
done
wait

#check if files exists
ls Processed_data/*/tophat/accepted_hits.bam &>/dev/null || ( echo "accepted_hits.bam file not found" && exit 2 )
echo 'tophat is done'

#sorting files in preparation for hts_count and hts_count
for i in Processed_data/*
do
si=$(echo $i | cut -d "/" -f2)

echo "Starting the samtools on ...." "$si"
samtools sort -on Processed_data/"$si"/tophat/accepted_hits.bam - | samtools view -h - | \
htseq-count -s no -m intersection-nonempty - \
/nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
> Processed_data/"$si"/"$si"_hts.txt

sample=$(echo 'sample_'$si)
sed "1s/.*/gene $sample/" "$i"/"$si"_hts.txt > Processed_data/"$si"/"$si"_htscounts.txt #file with named sample

echo "$si" "hts is done."
done
echo 'htscounts are complete!'

#part5. joining htscount output files into a single file in preparation for R processing
rm -f hts_count_merged.txt

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
[ -e "$i"/"$si"_htscounts.txt ] && awk '{print $1}' "$i"/"$si"_htscounts.txt > hts_count_merged.txt
done 
echo "Merge files were created for joining"

for i in Processed_data/* ; do
si=$(echo $i | cut -d "/" -f2)
#head -2 "$i"/"$si"_htscounts.txt
[ -e "$i"/"$si"_htscounts.txt ] && join hts_count_merged.txt "$i"/"$si"_htscounts.txt >tmp_merged.txt
echo "$si" "was hts joined"

mv tmp_merged.txt hts_count_merged.txt
done

echo "The script was finished on:" $(date)

