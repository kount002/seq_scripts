# data files obtained form Micro MiSeq run
# Sript used to cluster quality filtered and trimmed sequences
# will be used to calculate diversity
# This scrip is based on  phm_deversity/phm_diverity.sh that cleans the MiSeq reads and assembles them.


# Requires running phm_diversity.sh before
# START SCRIPT from within the directory
############################# 

# #part1. preprocessing QC: fastqc

#sed -n 'N,$p' phm_diversity.sh | bash

PDIR=$(pwd)


for i in Processed_data/*/tophat/accepted_hits.bam ; do
si=$(echo $i | cut -d "/" -f2i)
dirn=$(dirname "$i")
file=$(basename "$i")
clrt=$(dirname "$dirn")

mkdir -p "$clrt"/clustering


#convert accepted reads after tophat assembly to fasta file. Does not distinguish between R1 and R2
samtools view $i | awk '{OFS="\t"; print ">"$1"\n"$10}' - > $(dirn)/accepted.fasta

#pacard tool estimete lib size
java -jar /nfs/central/home/josh/collabs/software/picard-tools-1.135/picard.jar EstimateLibraryComplexity INPUT=accepted_hits.bam OUTPUT=pictest.txt



(
#sort seq by the length; required for further processing
uclust --sort $(dirn)/accepted.fasta  --output "$clrt"/clustering/"$si"_clean.sort.fasta


uclust --input "$clrt"/clustering/"$si"_clean.sort.fasta --uc "$clrt"/clustering/"$si"_clust.uc --id 0.95
uclust --uc2fasta "$clrt"/clustring/"$si"_clust.uc --input "$clrt"/clustering/"$si"_clean.sort.fasta \
--output "$clrt"/clustering/"$si"_out.fasta --types S
)&
done
wait
echo Cleanup complete!




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
mv tmp_merged.txt hts_count_merged.txt

done


