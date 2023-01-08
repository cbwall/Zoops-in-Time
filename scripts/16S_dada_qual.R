######### First in a series of scripts for R
# before you run make sure you close conda 
# in terminal: conda info --env 
# close open conda with 'conda deactivate'

# load dada2 package
library("devtools")
library("dada2")

start_time <- Sys.time() # track timing


#########

setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

# format metadata
run.metaD<- read.csv("/projects/ps-shurinlab/users/cbwall/Zoops_MP/data/Zoops_MP_metadata.csv")
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# ID the - controls
run.metaD$sample_control<-as.factor(run.metaD$Sample_Type)
run.metaD$sample_control<-gsub(",", "", run.metaD$sample_control) # delete commas in sample names
run.metaD$sample_control[run.metaD$sample_control=='Control PCR_NEG_1' |
                           run.metaD$sample_control=='Control PCR_NEG_2' |
                           run.metaD$sample_control=='Control PCR_NEG_3' |
                           run.metaD$sample_control=='Control PCR_NEG_4' |
                           run.metaD$sample_control=='Control PCR_NEG_5' |
                           run.metaD$sample_control=='Control PCR_NEG_6' |
                           run.metaD$sample_control=='Control PCR_NEG_7' |
                           run.metaD$sample_control=='Control PCR_NEG_8' |
                           run.metaD$sample_control=='Control DNA_CTRL_2' | 
                           run.metaD$sample_control=='Control DNA_CTRL_3' |
                           run.metaD$sample_control=='Control DNA_CTRL_4' |
                           run.metaD$sample_control=='Control DNA_CTRL_5' |
                           run.metaD$sample_control=='Control DNA_CTRL_6' | 
                           run.metaD$sample_control=='Control DNA_CTRL_7' |
                           run.metaD$sample_control=='Control WATER_DNA_CTRL_1' |
                           run.metaD$sample_control=='Control WATER_DNA_CTRL_2' |
                           run.metaD$sample_control=='DNA_CTRL_1'] <- "neg.controls"

run.metaD$sample_control[run.metaD$sample_control=='Control PCR_POS_1' |
                           run.metaD$sample_control=='Control PCR_POS_2' |
                           run.metaD$sample_control=='Control PCR_POS_3' |
                           run.metaD$sample_control=='Control PCR_POS_4'] <- "pos.controls"


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

write.csv(run.metaD, "output/run.metaD.edit.csv")

################################ Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#fnFs[1:3]

########################################
# create pdf of quality profiles for forward samples
# the aggregate = TRUE gives quality plot for ALL reads

pdf("output/qual_profiles_F.pdf")
plotQualityProfile(fnFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("output/qual_profiles_R.pdf")
plotQualityProfile(fnRs, aggregate = TRUE)
dev.off()

# final output is two pdfs of quality score profiles - forward and reverse
#################################

Sys.time() - start_time

