for i in 2 3 4
do
file=$(head -1 hts_count_hiseq_nano_merged.txt| cut -d " " -f $i)
tail -n +2 hts_count_hiseq_nano_merged.txt |\
head -n -5 hts_count_hiseq_nano_merged.txt |\
awk -v x=$i '$x>1000 {print $1}' > gene_counts_venn/"$file"_1000.txt
done
