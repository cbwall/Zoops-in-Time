

############## In terminal -- you need to do the following to get started in R

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
pacman::p_load('BiocManager', 'knitr', 'tidyverse', "dada2", "phyloseq", "decontam")

#############################################
##########################################################################################

# R code here for TSCC .sh job

# load dada2 package
library(dada2)
packageVersion('dada2')


start_time <- Sys.time() # track timing


#########
# format metadata
run.metaD<- read.csv("data/Zoops_MP_metadata.csv")
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# load in the cut-adapt samples in the "trimmed" folder
# run this in terminal 'gunzip data/trimmed_sequences/*'
# now unzipped... proceed

# path to folder containing demultiplexed library sequencing files
# make sure unzipped
miseq_path<-"data/trimmed" 
list.files(miseq_path)


################################# Filter and Trim
### remove low quality reads, trim to consistent length

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_trimmed.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_trimmed.fastq"))

sampleNames.p2 <- sapply(strsplit(fnFs, "_"), `[`, 2) # extract sample names
sampleNames.p3 <- sapply(strsplit(fnFs, "_"), `[`, 3) # extract the run # sample
sampleNames<-paste(sampleNames.p2,sampleNames.p3) # compile
sampleNames<-gsub(" ", "_", sampleNames) # remove space and add an underscore

# remove the "Undetermined" sample name 'S0_R1'
sampleNames<-sampleNames[1:368]


#### add this SampleNames to the metadata file
run.metaD$sampleNames<-sampleNames

################################ Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#fnFs[1:3]

############################## quality score plot for forward reads
# plotQualityProfile(fnFs[c(1,12)])

# quality score plot for reverse reads
# plotQualityProfile(fnRs[c(2,8)])


setwd("/projects/ps-shurinlab/cbwall/Zoops_MP/output")

# create pdf of quality profiles for forward samples
# tried to make a loop for this but it didn't work so this is hardcoded for 105 samples
# edit lines as needed for your number of samples - four at a time looks pretty good but can definitely be condensed
pdf("qual_profiles_F.pdf")
plotQualityProfile(fastqFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("qual_profiles_R.pdf")
plotQualityProfile(fastqRs, aggregate = TRUE)
dev.off()

Sys.time() - start_time

# final output is two pdfs of quality score profiles - forward and reverse


filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sampleNames, "_F_trimfilt.fastq"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_trimfilt.fastq"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,100), #trimLeft=c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

# this dataframe will show us the files and the # of reads going in, and out (post filterAndTrim)
head(out)

#inspect quality score plots for the filtered reads
plotQualityProfile(filtFs)
plotQualityProfile(filtRs)

# write to the export folder
write.csv(out, file="output/out.trim.csv")

############################ Let's learn and inspect the error rate
### estimate the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
# 171551200 total bases in 857756 reads from 33 samples will be used for learning the error rates.
saveRDS(errF, file="output/errF.rds")
errF<- readRDS("output/errF.rds")

errR <- learnErrors(filtRs, multithread=TRUE)
# 100007100 total bases in 1000071 reads from 80 samples will be used for learning the error rates.
saveRDS(errR, file="output/errR.rds")
errR<- readRDS("output/errR.rds")


# plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
