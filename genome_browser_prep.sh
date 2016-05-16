#check if files exists
ls Processed_data/*/tophat/accepted_hits.bam &>/dev/null || ( echo "accepted_hits.bam file not found" && exit 2 )
echo 'accepted file is present'

mkdir -p genbrowser_files

#sorting files in preparation for hts_count and hts_count
for i in Processed_data/*
do
si=$(echo $i | cut -d "/" -f2) #sample name

echo "Starting the samtools on ...." "$si"
samtools sort Processed_data/"$si"/tophat/accepted_hits.bam genbrowser_files/"$si"_accepted_sorted

samtools index genbrowser_files/"$si"_accepted_sorted.bam genbrowser_files/"$si"_accepted_sorted.bam.bai

echo "$si" "indexing is done."
done
