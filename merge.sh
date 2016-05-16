echo "$0" "script was started..." "$(date)"
rm -f out.merge.bam
i=~/anal/human/hiseq_trio/Processed_data/560/tophat/accepted_hits.bam
i2=~/anal/human/hiseq_trio/Processed_data/561/tophat/accepted_hits.bam
i3=~/anal/human/hiseq_trio/Processed_data/570/tophat/accepted_hits.bam

samtools merge out.merge.bam $i $i2 $i3
wait
echo "start preseq script"
bash run.preseq.sh out.merge.bam RNAMerge &> merge.log
wait
rm out.merge.bam
echo "Sript is finished" "$(date)"
