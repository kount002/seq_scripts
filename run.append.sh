#append files from two different runs
PDIR=$(pwd)

mkdir -p "$PDIR"/Processed_data/tmp

for i in ~/Seq_data/Human/Micro/Phgm_diversity/*R*.fastq.gz ; do

si=$(echo $i | basename $i | cut -c 1-5 )

[ -e  "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq ] && rm  "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq
[ -e  "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq ] && rm  "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq

(for j in ~/Seq_data/Human/Micro/Phgm_diversity*/$si*R1*.fastq.gz ; do
zcat -c $j >> "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq ; done
gzip "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq
rm -f "$PDIR"/Processed_data/tmp/"$si"_R1_tmp.fastq)&


(for j in ~/Seq_data/Human/Micro/Phgm_diversity*/$si*R2*.fastq.gz ; do
zcat -c $j >> "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq ; done
gzip "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq
rm -f "$PDIR"/Processed_data/tmp/"$si"_R2_tmp.fastq)&

done
wait
echo 'Raw sequence files are concatenated'

