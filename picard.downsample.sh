#runs picard downsamplesam command
# also tips on running samtool downsampling

if [ $1 ]
then
echo "Sample $1 of the reads"
fi
echo "Enter fraction to sample 0.0 - 1.0:"
read fraction
echo "Using $fraction"

java -jar /nfs/central/home/josh/collabs/software/picard-tools-1.135/picard.jar DownsampleSam PROBABILITY=$fraction INPUT=accepted_hits.bam OUTPUT=accepted_hits.downsample"$fraction".bam
