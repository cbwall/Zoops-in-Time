## download data from server

# TSCC login
ssh cbwall@tscc-login.sdsc.edu

# <enter password>

# change wd into the Shurin Lab projects and see the files

cd /projects/ps-shurinlab
pwd # show directory
ls # show files

### make new directory, navigate to the directory where you want new file to be
# this will house the raw data on the server
mkdir Yos_plank16S_time_MP
mkdir Horse_fish16S_MDJP

# made a mistake? rename file names? if empty...
# mv old_name new_name

# change directory to where you want MiSeq files to go
cd Yos_plank16S_time_MP

#######
#use wget to add to file location
#######
# Maddy 1 
wget -r -nH -nc -np -R index.html* "http://slimsdata.genomecenter.ucdavis.edu/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned1/Project_JSJD_MP_Pool1/"

# Maddy 2
wget -r -nH -nc -np -R index.html* "http://slimsdata.genomecenter.ucdavis.edu/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned2/Project_JSJD_MP_Pool1/"
#######


# change directory to where you want MiSeq files to go
cd /projects/ps-shurinlab/Horse_fish16S_MDJP

#######
# Josh 1
wget -r -nH -nc -np -R index.html* "http://slimsdata.genomecenter.ucdavis.edu/Data/eb80e4guok/221216_M01533_0053_MS3148640-500V2/Unaligned1/Project_JSJD_M1300P_Dominguez/"

# Josh 2
wget -r -nH -nc -np -R index.html* "http://slimsdata.genomecenter.ucdavis.edu/Data/eb80e4guok/221216_M01533_0053_MS3148640-500V2/Unaligned2/Project_JSJD_M1300P_Dominguez/"
#######


############## For Maddy,
# imported data lives in these 2 directories
# /Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned1/Project_JSJD_MP_Pool1
# /Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned2/Project_JSJD_MP_Pool1

# copy the files above into a single folder in user folder
# can copy (cp) files from one place to another with directory (-r)
cp -r /projects/ps-shurinlab/Trophobiomes/ANL.sequencing raw.sequences/

# copy from a folder based on a string into a new folder
cp CBW_YoTB_6[0-9][0-9]* For.Josh/

# navigate to user directory
cd /projects/ps-shurinlab/users/cbwall

# make new folder for the project
mkdir Zoops_MP

# make a folder for sequences
mkdir Zoops_MP/data
mkdir Zoops_MP/data/raw_sequences_MP

# show files in the main projects (first set of files, there are 2)
ls Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned1/Project_JSJD_MP_Pool1
ls Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned2/Project_JSJD_MP_Pool1

# move files ( the /* says move contents of the files in the "Project_JSJD_MP_Pool1" folder)
# from old dir to new dir
mv -v /projects/ps-shurinlab/Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned1/Project_JSJD_MP_Pool1/* /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/raw_sequences_MP/
mv -v /projects/ps-shurinlab/Yos_plank16S_time_MP/Data/mxhmt2vqki/221216_M02034_0667_MS3148658-500V2/Unaligned2/Project_JSJD_MP_Pool1/* /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/raw_sequences_MP/

# navigate to WD and remove any index or other files
cd Zoops_MP/data/raw_sequences_MP/ 
rm laneBarcode.html
rm @md5Sum.md5 

# add in metadata csv file to the 'Zoops_MP' folder
# must open a new terminal window on your home directoy, will ask for password
scp ~/Downloads/Zoops_MP_metadata.csv cbwall@tscc-login.sdsc.edu:/projects/ps-shurinlab/users/cbwall/Zoops_MP

# make folders and subfolders for trimmed and filtered data in the "Zoops_MP" folder
mkdir data/trimmed
mkdir output

mkdir data/filtered
mkdir output
mkdir output/errs
mkdir output/outs

mkdir output/outs/cutadapt_out
mkdir output/errs/cutadapt_err

# create text file with list of all sample names
# assumes samples are in folder data with month.year after, and have format of...
# sampleid_S###_L001_R1_001.fastq
ls /projects/ps-shurinlab/users/cbwall/Zoops_MP/data/raw_sequences_MP | sed 's/_L[[:digit:]]\+_.*//g' > samples.txt

# see if it worked
head samples.txt
cat samples.txt # see whole file

################# run the script in TSCC
# move shell (.sh) script from your machine to TSCC for job processing
scp ~/Downloads/cutadapt_16S_MPzoops.sh cbwall@tscc-login.sdsc.edu:/projects/ps-shurinlab/users/cbwall/scripts


#/#/#/#
# ready to submit job

# used this for MP
# navigate to user folder where 'scripts' is then run...

qsub scripts/cutadapt_16S_MPzoops.sh
# will give a job ID, like '30508501.tscc-mgr7'
# to check status, run 
qstat  30524649

#delete a job
qdel <jobid>

# check job status
checkjob <jobid>

# check nodes for a user
qstat -n -u cbwall

# inspect the output file to confirm no error reports

#/#/#/##/#/#/##/#/#/##/#/#/##/#/#/#
 # success! Now confirm it
 
cd data/raw_sequences_MP/ # change WD to see original, raw files
gunzip *.fastq.gz # unzip

# pick some files at random to inspect
# this shows  occurrences in the orignal
grep "GTGCCAGCCGCCGCGGTAA" MP_171_16S_S171_L001_R1_001.fastq | wc -l  # gives 56
grep "GTGCCAGCCGCCGCGGTAA" MP_100_16S_S100_L001_R1_001.fastq | wc -l # gives 8
grep "GTGCCAGCCGCCGCGGTAA" MP_08_16S_S8_L001_R1_001.fastq | wc -l # gives 19

# rezip 
gzip *.fastq

# check the trimmed files
cd data/trimmed/ # change WD to see original, raw files
gunzip *.fastq.gz # unzip

grep "GTGCCAGCCGCCGCGGTAA" MP_171_16S_S171_R1_trimmed.fastq | wc -l  # gives 0
grep "GTGCCAGCCGCCGCGGTAA" MP_100_16S_S100_R1_trimmed.fastq | wc -l # gives 0
grep "GTGCCAGCCGCCGCGGTAA" MP_08_16S_S8_R1_trimmed.fastq | wc -l # gives 0

# SUCCESS!
##################

# export the any files to view or inspect (or run if you can manage) on home device. do this from your home location in terminal
scp -r cbwall@tscc-login.sdsc.edu:/projects/ps-shurinlab/users/cbwall/Zoops_MP/data/trimmed/*  ~/Downloads/Zoops_MP



###### moving on to R in TSCC .....


############## R
module load R 
R --version

# to get new version (4.1.2)
export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.1.2
R --version

# install and load packages
# will ask you to install into a personal library since this one not writable, say yes x2

export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.1.2
R # now R is loaded

install.packages('pacman')

# use pacman to load CRAN packages missing
pacman::p_load('knitr', 'tidyverse', "dada2", "phyloseq", "decontam")

# to exit R module use 
# quit()


# All necessary packages should be installed before you submit a job that runs an R script.
# To submit a job that runs the script 'mycode.R', create and save a .sh file:

# it would be a good idea to have your R project directory being similar in structure and naming as your TSCC folders
# you are now ready to run an R script in TSCC


