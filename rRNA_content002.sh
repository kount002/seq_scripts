# modified from phm_diversity_genelist.sh, 
#This script cleanups bacterial sequences and run them against SortMeRNA to evaluate the proportion of rRNA in the samples prepared at the core with use of Ribozero.

#The script will output log file and two files with rRNA sequences and non_rRNA
 
#usage: Script.sh path_to_directory_with_one_sample


############################# 

#part 1: verify agruments presence

if [ -z $1 ] ; then
echo "Usage: rRNA_content.sh path_to_directory_with_one_sample"
exit 2
fi


#part 2. prepare the environment for the processing.
echo "The script_" $(basename $0) "_was started at:" $(date)

PDIR=$(pwd)

#make processing directory to place work files
mkdir -p /nfs/gems_sata/tedder/evgueni/anal/tmp_rRNA
[ -e Processed_data ] || ln -s /nfs/gems_sata/tedder/evgueni/anal/tmp_rRNA/ Processed_data


#part 3. concatenate files
si=$(basename "$1")
dir=$(echo $si | cut -c 1-4)
echo "Processing" $dir

rm -f Processed_data/"$dir"/*
mkdir -p Processed_data/"$dir"
touch Processed_data/"$dir"/"$dir"_R1.fastq
touch Processed_data/"$dir"/"$dir"_R2.fastq
(
for i in "$1"/*R1* ; do
echo $i "test"
zcat $i >> Processed_data/"$dir"/"$dir"_R1.fastq
done) &

(
for i in "$1"/*R2* ; do
zcat $i >> Processed_data/"$dir"/"$dir"_R2.fastq
done) &
wait
echo "Done concatenate" $(ls Processed_data/"$dir")


#compress files
for i in Processed_data/"$dir"/* ; do
(gzip $i)&
done
wait

echo "Done with compression"


#part 4. adapter removal: cutadapt alone or wt combination with Trim Galore
#run cutadapt
file1=$(ls Processed_data/"$dir"/*R1*)
file2=$(ls Processed_data/"$dir"/*R2*)

# R1 reads use universal primers for sequencing
#for R1 look for complement of 3-end of insert + Index adapter: GTTGCGGCCGCTGGATTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ( insert sequence starts with CCATGGCCGCCGAGAAC edd to universal adapter in reverse when cleaning for R2)
#Need to remove second adapter at 5' end that is insert start seq
~/.local/bin/cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz -p "$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz \
$file1 $file2 \
&> Processed_data/"$dir"/"$dir"_cutadapt.log

#for R2 look for reverse complement of insert start and Universal primer (GTTCTCGGCGGCCATGG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT)
~/.local/bin/cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q20 -m20 -e 0.1 \
-o "$PDIR"/Processed_data/"$dir"/"$dir"_R2_trimmed.fastq.gz -p "$PDIR"/Processed_data/"$dir"/"$dir"_R1_trimmed.fastq.gz \
"$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz \
&>> Processed_data/"$dir"/"$dir"_cutadapt.log

rm "$PDIR"/Processed_data/"$dir"/"$dir"_atmp1.fastq.gz "$PDIR"/Processed_data/"$dir"/"$dir"_atmp2.fastq.gz

wait
echo 'All adapters are removed'


#decompress files
for i in Processed_data/"$dir"/* ; do
(gunzip $i)&
done
wait

#interleaving paired files into single read file
bash /winhomes/kount002/binaries/sortmerna-2.0-linux-64/scripts/merge_paired_reads.sh Processed_data/"$dir"/"$dir"_R1_trimmed.fastq  Processed_data/"$dir"/"$dir"_R2_trimmed.fastq  Processed_data/"$dir"/"$dir"_interleaved.fastq
echo "Done with uncompression"

#part4. Run SortmeRNA using supplied library

/winhomes/kount002/binaries/sortmerna-2.0-linux-64/sortmerna --ref \
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-16s-id90.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-bac-16s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-23s-id98.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-bac-23s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-16s-id95.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-arc-16s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-23s-id98.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-arc-23s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-18s-id95.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-euk-18s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-28s-id98.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/silva-euk-28s:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/rfam-5s-database-id98.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/rfam-5s-db:\
/winhomes/kount002/binaries/sortmerna-2.0-linux-64/rRNA_databases/rfam-5.8s-database-id98.fasta,/winhomes/kount002/binaries/sortmerna-2.0-linux-64/index/rfam-5.8s-db \
--reads Processed_data/"$dir"/"$dir"_interleaved.fastq --paired_in --fastx --aligned Processed_data/"$dir"/"$dir"_rRNA --other Processed_data/"$dir"_non_rRNA --log -a 6 -m 4096 -v


bash /winhomes/kount002/binaries/sortmerna-2.0-linux-64/scripts/unmerge_paired_reads.sh Processed_data/"$dir"/"$dir"_non_rRNA.fastq \
Processed_data/"$dir"/"$dir"_R1_non_rRNA.fastq  Processed_data/"$dir"/"$dir"_R2_non_rRNA.fastq

echo "Sorting is done"

echo "Script is finished " $(date)

