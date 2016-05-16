# the script contatenates files from the poorly binned Illumina ran and then attempts to demultiplex them
#Usage 
#Map file contains tabulated Index and sample_name

#check if file was provided with a argument 
if [ -z "$1" ] || [ -z "$2" ]; then
echo "Sript usage: " "$0" "path_dir_with_raw_data  " "map_file"
exit
fi

echo "$0" "script has been started" "$(date)"

sample=$(basename $1)

#Part 1 append files from tree runs to one big file in the processing directory
mkdir -p "$1"tmp
#remove preexisting files
rm -f "$1"tmp/*

#append cycles for R1 and R2
(for j in "$1"*_R1_* ; do
zcat $j >> "$1"tmp/all_R1_app.fastq ; done)&
echo "R1 is processing"

(for j in "$1"*_R2_* ; do
zcat $j >> "$1"tmp/all_R2_app.fastq ; done)&
echo "R2 is processing"

(for j in "$1"*_I1_* ; do
zcat $j >> "$1"tmp/all_I1_app.fastq ; done)&
echo "I1 is processing"

wait
echo 'Raw sequence files are concatenated'

mkdir -p "$1"demultiplexed

(fastq-multx -g "$1"tmp/all_I1_app.fastq "$1"tmp/all_R1_app.fastq "$1"tmp/all_R2_app.fastq -o "$1"demultiplexed/%_R1.fastq -o "$1"demultiplexed/%_R2.fastq) && rm -r "$1"tmp

#(iu-demultiplex -s "$1""$2" --r1 "$1"tmp/all_R1_app.fastq --r2 "$1"tmp/all_R2_app.fastq --index all_I1_app.fastq -o "$1"demultiplexed) && rm -r "$1"tmp
echo "Finished demultiplexing"


echo "The script was finished.." "$(date)"


