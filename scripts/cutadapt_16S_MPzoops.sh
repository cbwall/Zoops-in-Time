#!/bin/bash
#PBS -q hotel
#PBS -N cutadapt
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:30:00
#PBS -o /projects/ps-shurinlab/users/cbwall/Zoops_MP/output/outs/cutadapt_out
#PBS -e /projects/ps-shurinlab/users/cbwall/Zoops_MP/output/errs/cutadapt_err
#PBS -V
#PBS -M cbwall@ucsd.edu
#PBS -m ae
#PBS -A shurin-group


# activate conda and cutadapt, if there are problems the culprit is probably one of these two lines
source ~/.bashrc
conda activate cutadaptenv


# change directory to easily access fastq sequence files
cd /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/raw_sequences_MP

# create text file with list of all sample names
# assumes samples are in folder data with month.year after, and have format of...
# sampleid_S###_L001_R1_001.fastq
ls /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/raw_sequences_MP | sed 's/_L[[:digit:]]\+_.*//g' > /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/samples.txt


#remove any index files
rm index*
rm laneBarcode.html
rm @md5Sum.md5 


# unzip sequence files, if necessary
# gunzip *.fastq.gz


# change to lab directory and create output file in personal folder
# cd /projects/ps-shurinlab/users/cbwall/Zoops_MP

# --- back to the code...
# code written specifically for samples that only have the 16S adapter sequence on the 5' end 
# (this should be standard for our library runs)

### re: -g and -G
### this is 'regular 5' adapter' that is not anchored
### 5’ adapter is a piece of DNA ligated to the 5’ end of the DNA fragment of interest.
### For this type of adapter to be found, the adapter sequence needs to either appear in full
### somewhere within the read (internal match) or at the start (5’ end) of it
### In all cases, the adapter itself and the sequence preceding it is removed.

# -g ADAPTER for forward primer/adapter
# -G adapter for reverse primer/adapter
# -o for output of forward read
# -p for output of reverse read
# -m for minimum and -M for maximum read length

# works for samples with the file names in format: SAMPLEID_SAMPLENUMBER_L001_R(1/2)_001.fastq, 
# may need to adjust

# also, cutadapt output stats (“cutadapt_primer_trimming_stats.txt”), assess how things went

#### --- back to the code...


#run the loop


for SAMPLEID in $(cat data/samples.txt);
do
    echo "On sample: $SAMPLEID"
    cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT \
    -o data/trimmed/${SAMPLEID}_R1_trimmed.fastq \
    -p data/trimmed/${SAMPLEID}_R2_trimmed.fastq \
    data/raw_sequences_MP/${SAMPLEID}_L001_R1_001.fastq \
    data/raw_sequences_MP/${SAMPLEID}_L001_R2_001.fastq \
	>> output/cutadapt_primer_trimming_stats.txt 2>&1
done



# re-zip sequence files
gzip data/raw_sequences_MP/*.fastq
gzip data/trimmed/*.fastq

# final output is a set of trimmed fastq files in a "trimmed" folder
