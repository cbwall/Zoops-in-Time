############ 16S part 2

# load dada2
library(dada2)
packageVersion('dada2')

#start_time <- Sys.time() # track timing


###################################################################################################
#### need the below if planning to run this script solo, otherwise combine in a single shell
#################################
#setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

# read in formatted metaData
run.metaD<-read.csv("output/dada proc/run.metaD.edit.csv")
run.metaD<- run.metaD[-1]
row.names(run.metaD)<-run.metaD$sampleNames

# Sample Names
sampleNames<-run.metaD$sampleNames

# Sort ensures forward/reverse reads are in same order
miseq_path<-"data/trimmed" 
list.files(miseq_path)

fnFs <- sort(list.files(miseq_path, pattern="_R1_trimmed.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_trimmed.fastq"))

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

filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(225,200), #trimLeft=c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

# this dataframe will show us the files and the # of reads going in, and out (post filterAndTrim)
# head(filt_out)

#################################
# write to the export folder
#write.csv(filt_out, file="/projects/ps-shurinlab/users/cbwall/Zoops_MP/output/filt_out.csv")
write.csv(filt_out, file="output/dada proc/filt_out.csv")
filt_out<-read.csv("output/dada proc/filt_out.csv")
row.names(filt_out)<-filt_out$X
filt_out<-filt_out[-1]

#################################
#inspect quality score plots for the filtered reads
pdf("output/dada proc/qual_profiles_filtF.pdf")
plotQualityProfile(filtFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("output/dada proc/qual_profiles_filtR.pdf")
plotQualityProfile(filtRs, aggregate = TRUE)
dev.off()


########################## learn and inspect the error rate
# estimate the error rates F
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF, "output/dada proc/errF.rds")
errF<- readRDS("output/dada proc/errF.rds")

# error rate for reverse
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errR, "output/dada proc/errR.rds")
errR<- readRDS("output/dada proc/errR.rds")

# plot error rates
pdf("output/dada proc/F_err.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("output/dada proc/R_err.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()


########################## Dereplicate 
### Derep
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "output/dada proc/derepFs.rds")
#derepFs<- readRDS("output/TSCC-output/derepFs.rds")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs, "output/dada proc/derepRs.rds")
#derepRs<- readRDS("output/TSCC-output/derepRs.rds")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

########################## The DADA2 sequence inference method can run in two different modes: 
# By default, the dada function processes each sample independently. However, pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples. The dada2 package offers two types of pooling. dada(..., pool=TRUE) performs standard pooled processing, in which all samples are pooled together for sample inference. dada(..., pool="pseudo") performs pseudo-pooling, in which samples are processed independently after sharing information between samples, approximating pooled sample inference in linear time.

# pool= TRUE and pseudo, we will use TRUE here
# TRUE takes a lot more time and memory
dadaFs <-dada(derepFs, err=errF, multithread=2, pool="pseudo")
saveRDS(dadaFs, "output/dada proc/dadaFs.rds")
#dadaFs<- readRDS("output/TSCC-output/dadaFs.rds")

dadaRs <-dada(derepRs, err=errF, multithread=2, pool="pseudo")
saveRDS(dadaRs, "output/dada proc/dadaRs.rds") 
#dadaRs <-readRDS("output/TSCC-output/dadaRs.rds") 

# inspect data
dadaFs[[1]]

########################## merge the pairs of dadaFs (ASVs F), dereplicated-F, dadaRs (ASVs R), and dereplicated-R
# trim overhang in case any reads go past opposing primers
# can also set 'minOverlap=' if there is concern for low overlap -- not the case here with the PE-300
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang=TRUE, verbose=TRUE)
saveRDS(mergers, "output/dada proc/merged_amplicons.rds")
#mergers<-readRDS("output/TSCC-output/merged_amplicons.rds")

# remove anything that is "mock community"
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "output/dada proc/seqtab.rds")
#seqtab <- readRDS("output/TSCC-output/seqtab.rds")

####
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants
pdf(file="output/dada proc/seqtab.hist.pdf", height=5, width=5)
hist(nchar(getSequences(seqtab)))
dev.off()

table(nchar(getSequences(seqtab)))
dim(seqtab) #368 by 29278

### remove chimera
seqtab.nochim <-removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "output/dada proc/seqtab.nochim.rds")
#seqtab.nochim<-readRDS("output/TSCC-output/seqtab.nochim.rds")

# how many samples cleared processing
sum(seqtab.nochim)/sum(seqtab) # 96% of samples kept, chimeras ~ 3% of merged reads


###################### summary table
# create table showing read count loss throughout the process
getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(row.names = sampleNames,
                          dada2_input = filt_out[,1], filtered = filt_out[,2],
                          dada_f = sapply(dadaFs, getN), dada_r = sapply(dadaRs, getN),
                          merged = sapply(mergers, getN),
                          nonchim = rowSums(seqtab.nochim),
                          final_perc_reads_retained = round(rowSums(seqtab.nochim)/filt_out[,1]*100,1))

write.table(summary_tab, "output/dada proc/read_count_tracking_local.tsv", quote = FALSE, sep = "\t", col.names = NA)

##### assign taxonomy using SILVA 2021 database
taxa <- assignTaxonomy(seqtab.nochim, "output/SILVA/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# to add in Species for 16S
taxa <- addSpecies(taxa, "output/SILVA/silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "output/taxa.rds")
#taxa<-readRDS("output/TSCC-output/taxa.rds")
  
# final output is folder of filtered reads, summary table, 10 .rds files

