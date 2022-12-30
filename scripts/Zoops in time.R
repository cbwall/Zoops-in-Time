


######## load packages

if (!require('knitr')) install.packages('knitr'); library('knitr')
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center')

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load CRAN packages missing
pacman::p_load('knitr', 'tidyverse', 'knitr', 'magrittr', 'effects', 'devtools',
               'stringi', 'dplyr', "ggplot2", "gridExtra", "dada2", "phyloseq", "vegan", "cowplot",
               "decontam","BiocManager", "dada2", "decipher")


#########
# format metadata
run.metaD<- read.csv("data/Zoops_MP_metadata.csv")
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# load in the cut-adapt samples in the "trimmed" folder
# run this in terminal 'gunzip data/trimmed_sequences/*'
# now unzipped... proceed

miseq_path<-"data/trimmed" # containing the fastq files after unzipping.
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
fnFs[1:3]

############################## quality score plot for forward reads
plotQualityProfile(fnFs[c(1,12)])

# quality score plot for reverse reads
plotQualityProfile(fnRs[c(2,8)])

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

########################## Dereplicate 
### Derep
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, file="output/derepFs.rds")
derepFs<- readRDS("output/derepFs.rds")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs, file="output/derepRs.rds")
derepRs<- readRDS("output/derepRs.rds")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
 
############################### #The DADA2 sequence inference method can run in two different modes: pool= TRUE and pseudo, we will use TRUE here
#####
dadaFs <-dada(derepFs, err=errF, multithread=2, pool=TRUE)
saveRDS(dadaFs, file="output/dadaFs.rds")
# shows 368 samples were pooled: 8,475,897 reads in 172,0334 unique sequences.

dadaRs <-dada(derepRs, err=errF, multithread=2, pool=TRUE)
saveRDS(dadaRs, file="output/dadaRs.rds")
# 

# inspect data
dadaFs[[1]]

############### merge the pairs of dadaFs (ASVs F), dereplicated-F, dadaRs (ASVs R), and dereplicated-R
# trim overhang in case any reads go past opposing primers
# can also set 'minOverlap=' if there is concern for low overlap -- not the case here with the PE-300
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang=TRUE, verbose=TRUE)

# remove anything that is "mock community"
seqtab <- makeSequenceTable(mergers)

#################################### if wanted to remove "mock community" (pcr_pos) you could do this here... we will do it later so we can evaluate the mock community
# seqtab <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

table(nchar(getSequences(seqtab)))
dim(seqtab) # shows 2437 ASVs in 18 samples


################################### remove chimera

seqtab.nochim <-removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# infers 459 chimeras

# how many samples cleared processing
sum(seqtab.nochim)/sum(seqtab) # 95% of samples kept, chimeras ~ 3% of merged reads

###################### summary table
getN <-function(x)sum(getUniques(x))
summary_tab <-cbind(out,sapply(dadaFs, getN), 
                    sapply(dadaRs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(summary_tab) <-c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(summary_tab) <- sampleNames
head(summary_tab)

write.csv(summary_tab, "output/read-count-tracking.csv")


#################################### assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "data/silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# to add in Species for 16S
taxa <- addSpecies(taxa, "data/silva/silva_species_assignment_v138.1.fa.gz")

# inspect taxonomic assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# let's save this as .RData since it is so time consuming!
saveRDS(taxa, file="output/taxaTable.rds")



##### #If the `assignTaxonomy` with SILVA is a pain, let's just load in the R data and run it. This is what was generated hen I ran the code above.

taxa<- readRDS("output/taxaTable.rds")
# now you've loaded in the code chunk info above and can proceed