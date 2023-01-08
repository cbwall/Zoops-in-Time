############ 16S part 2

# load dada2
library(dada2)
packageVersion('dada2')

start_time <- Sys.time() # track timing


###################################################################################################
#### need the below if planning to run this script solo, otherwise combine in a single shell
#################################
setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

# format metadata
run.metaD<- read.csv("/projects/ps-shurinlab/users/cbwall/Zoops_MP/data/Zoops_MP_metadata.csv")
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# Sort ensures forward/reverse reads are in same order
miseq_path<-"data/trimmed" 
list.files(miseq_path)

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

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
###################################################################################################
# end required loading


### Filter and Trim
filt_path <- "data/filtered" # Place filtered files in filtered/ subdirectory
#if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sampleNames, "_F_trimfilt.fastq"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_trimfilt.fastq"))


################# learn and inspect the error rate
# estimate the error rates F
errF<- readRDS("output/errF.rds")

# error rate for reverse
errR<- readRDS("output/errR.rds")

########################## Dereplicate 
### Derep
derepFs<- readRDS("output/derepFs.rds")
derepRs<- readRDS("output/derepRs.rds")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

## The DADA2 sequence inference method can run in two different modes: pool= TRUE and pseudo, we will use TRUE here
#dadaFs <-dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaFs<- readRDS("output/dadaFs.rds")
# 368 samples were pooled: 8371792 reads in 1363837 unique sequences.

dadaRs <-readRDS("output/dadaRs.rds") 

# inspect data
dadaFs[[1]]
# 786 sequence variants were inferred from 3622 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

############### merge the pairs of dadaFs (ASVs F), dereplicated-F, dadaRs (ASVs R), and dereplicated-R
# trim overhang in case any reads go past opposing primers
# can also set 'minOverlap=' if there is concern for low overlap -- not the case here with the PE-300
mergers<-readRDS("output/merged_amplicons.rds")

# remove anything that is "mock community"
seqtab <- readRDS("output/seqtab.rds")

####
# inspect sequence length
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. This table contains 293 ASVs, and the lengths of our merged sequences all fall within the expected range for this V4 amplicon.

table(nchar(getSequences(seqtab)))
dim(seqtab) 
# 21,243 ASVs in 368 samples

### remove chimera
seqtab.nochim<-readRDS("output/seqtab.nochim.rds")

# how many samples cleared processing
sum(seqtab.nochim)/sum(seqtab) # 96% of samples kept, chimeras ~ 4% of merged reads


# read back in the filt_out
###################### summary table
filt_out<-read.csv("output/filt_out.csv")
# note below that it filt_out is loaded in (not pulling back in) then filt_out[,1] and [,2] works instead of 2,3
# create table showing read count loss throughout the process

getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(row.names = sampleNames,
                          dada2_input = filt_out[,2], filtered = filt_out[,3],
                          dada_f = sapply(dadaFs, getN), dada_r = sapply(dadaRs, getN),
                          merged = sapply(mergers, getN),
                          nonchim = rowSums(seqtab.nochim),
                          final_perc_reads_retained = round(rowSums(seqtab.nochim)/filt_out[,2]*100,1))
write.table(summary_tab, "output/read_count_tracking.tsv", quote = FALSE, sep = "\t", col.names = NA)

##### assign taxonomy using SILVA 2021 database
taxa <- assignTaxonomy(seqtab.nochim, "/projects/ps-shurinlab/Databases/SILVA/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# to add in Species for 16S
taxa <- addSpecies(taxa, "/projects/ps-shurinlab/Databases/SILVA/silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "output/taxa.rds")

# final output is folder of filtered reads, summary table, 10 .rds files

