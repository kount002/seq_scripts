join hts_count_merged.txt ../nanoset/hts_count_merged_nano.txt > tmp.txt
awk '{print $1, $2+$5, $3+$6, $4+$7}' tmp.txt > hts_count_hiseq_nano_merged.txt
rm tmp.txt
