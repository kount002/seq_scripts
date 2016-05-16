# same as below but for sequence set 1
#script to sort out phiX sequences to determine if the read1 file contains sequences from the custom primer among the phiX and stock Illumina primers

#the script is preceeded by command: 
#zcat Undetermined_S0_L001_R1_001.fastq.gz | fastq_quality_filter -Q33 -q 20 -p 90 > Undetermined.custom.filter.fastq
# can filter out sequences with stock illumina by grep -v '^CCATGG'
# after sequnces were cleaned with bowtie




bowtie2 -p 8 -x /nfs/gems_sata/tedder/evgueni/Reference/PhiX174/bowtie2_index/bt2_phix \
/home/kount002/Seq_data/Human/Micro/Phgm_diversity_2/"$1" > "$1"_output.sam --un "$1"_unmapped.fastq
