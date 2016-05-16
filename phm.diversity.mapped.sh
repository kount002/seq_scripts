# Estimates diversity using reads mapped to the human genome and preseq lc-extrapolation
# Builds on the phm_diversity_append.sh after it assembles all three samples from three first reads.



for i in Processed_data/*/tophat/accepted_hits.bam ; do
si=$(echo "$i" | cut -d "/" -f 2)
echo "Picked " "$i" "for samtools processing" "Use" "$si"
samtools sort "$i" Processed_data/tmp."$si".sorted
done

rm -f Processed_data/out.merged*.bam
samtools merge -f Processed_data/out.merged.bam Processed_data/tmp*sorted*
samtools sort Processed_data/out.merged.bam Processed_data/out.merged.sorted
echo "Finished with BAM file"
rm -f Processed_data/tmp*sorted*


#preseq requires files to be sorted

~/binaries/preseq-1.0.2.Linux_x86_64/preseq lc_extrap -B -o lc.extrap.appended.files.txt Processed_data/out.merged.sorted.bam

~/binaries/preseq-1.0.2.Linux_x86_64/preseq c_curve -B -v Processed_data/out.merged.sorted.bam &> c_curve.merged.txt

