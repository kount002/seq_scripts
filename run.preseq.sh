#scrip runs preseq lc-extrapolate and C-curve
#Usage
##	bash run.preseq.sh infile.bam tag 
#part 1: request input file 

##echo "Enter analysed BAM file"
##read file
#Usage run.preseq.sh path_to_accepted.reads.bam sample_label
if [ -z $1 ] || [ -z $2 ] ; then
echo "Usage run.preseq.sh path_to_accepted.reads.bam sample_label"
exit
fi

echo "$0" "sript was started" "$(date)"
file=$1
#perform sorting 
lfile=$(basename $file)
dir=$(dirname $file)

##############
# Test to run bed format conversion

samtools view -F 0x0004 $file |	awk '{OFS="\t"; if (and($2, 16)) \
		print $3,$4,$4+length($10),$1,$5,"-"; \
		else print $3,$4,$4+length($10),$1,$5,"+" }' > "$dir"/out_"$2".bed
echo "Finished BED conversion of ..." "$file"

#samtools view "$file" | \
 # awk '{split ($6,a,"[MIDNSHP]"); bp=$4-1; n=0; \
  #  for (i=1; i<=length(a); i++) { \
   #   n+=1+length(a[i]); \
    #  if (substr($6,n,1)=="M") print $3"\t"bp"\t"(bp+=a[i]); \
     # if (substr($6,n,1)=="D") bp+=a[i]; \
   # } \
 # }' > out.bed

sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 "$dir"/out_"$2".bed > "$dir"/out.sort_"$2".bed
rm "$dir"/out_"$2".bed
echo "Finished sorting of ..." "$file"

##############

#samtools sort "$file" tmp."$lfile".sorted
##java -Djava.io.tmpdir=$(pwd)/tmp -jar /nfs/central/home/josh/collabs/software/picard-tools-1.135/picard.jar SortSam \
##I="$file" O=tmp_"$lfile"_sorted.bam SO=coordinate TMP_DIR=$(pwd)/tmp

##echo "Finished with sorting" $(date)
#preseq requires files to be sorted; run preseq

~/binaries/preseq-1.0.2.Linux_x86_64/preseq lc_extrap -P -o "$dir"/lc.extrap."$2".txt "$dir"/out.sort_"$2".bed

~/binaries/preseq-1.0.2.Linux_x86_64/preseq c_curve -P -v "$dir"/out.sort_"$2".bed &> "$dir"/c_curve."$2".txt

echo "The script has finished at" "$(date)"
#rmdir tmp
#rm -f tmp*sorted*
