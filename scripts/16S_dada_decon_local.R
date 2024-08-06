############ 16S part 4
# decontam

# remember to run through the 4.1.2 version of R instead of default

#install.packages("ggplot2")
library("ggplot2")

#install.packages("decontam")
library("decontam")
packageVersion("decontam")

library("BiocManager")

#BiocManager::install("phyloseq")
library("phyloseq")
packageVersion("phyloseq")

#BiocManager::install("Biostrings", force=TRUE)
library("Biostrings")
packageVersion("Biostrings")

### the below help with formatting, editing, exporting phyloseq
#install.packages("remotes")
#remotes::install_github("adrientaudiere/MiscMetabar")
#remotes::install_github("peterolah001/BiMiCo")
library("MiscMetabar")
library("BiMiCo")


#start_time <- Sys.time() # track timing

#setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

## required load in
###################################################################################################
# read in formatted metaData
run.metaD<-read.csv("output/reanalysis/run.metaD.edit.csv")
run.metaD<-run.metaD[-1] # remove the "X" column

# format metadata (have to redo this if loading back in)
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake", "Plate", "Plate_name", "Well", "sample_control")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# Sample Names
sampleNames<-run.metaD$sampleNames

taxa <- readRDS("output/reanalysis/taxa.rds")
seqtab.nochim <- readRDS("output/reanalysis/seqtab.nochim.rds")

###################################################################################################
### end required load in 

# load the three files generated in tables.R
counts <- read.table(file = 'output/reanalysis/ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
tax_table <- read.table(file = 'output/reanalysis/ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
asv_fasta <- readRDS("output/reanalysis/ASV_fasta.rds")

# sanity check to make sure the tables loaded right
# counts colnames should be a list of the samples, tax_table colnames should be kingdom phylum class etc...
print("Counts columns")
colnames(counts)

# need this to remove "X" if reading back in
colnames(counts) <- gsub(x = colnames(counts), pattern = "X", replacement = "")
colnames(counts)

print("Tax columns")
colnames(tax_table)


## sample data
# metadata is run.metaD
all(rownames(seqtab.nochim) %in% run.metaD$sampleNames)

rownames(run.metaD) <- run.metaD$sampleNames

ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
              sample_data(run.metaD), 
              tax_table(taxa))
# 23978 taxa in 368 samples

# make a string of DNA names and add to phyloseq
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

###########################
# Show available ranks in the dataset
rank_names(ps)

table(tax_table(ps)[, "Phylum"], exclude = NULL)
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
# 348 Archaea, 21171 Bacteria, 191 Eukaryota 2268 NA

# remove NAs in taxonomic table
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # only Archaea and Bacteria
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Chloroplast"))

# just in case!, 9471 ASVs
ps<- subset_taxa(ps, Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class)) 

# chloroplasts sneaking in under Order in some cyanobacteria, now 8107 ASVs
ps <- subset_taxa(ps, Order!="Chloroplast" | is.na(Class)) 

# could start with removing ALL non bacterial sequences
# ps<-rm_nonbac(ps)
# by avoiding this step we are keeping in the Archaea

#### summary
# #8107 taxa in 368 samples

###########################
## ID contaminants
sample_data(ps)$is.neg <- sample_data(ps)$sample_control == "neg.controls" 
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) # which are contaminants? 8105 NO, 2 YES
head(which(contamdf.prev$contaminant))

###########################
### prune controls and low reads
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam
#8105 taxa in 368 samples

#make sure those negative and pos controls are out! 
ps.noncontam.controls.out<- subset_samples(ps.noncontam, 
                                           !(sample_control %in% "neg.controls"))
ps.noncontam.controls.out<- subset_samples(ps.noncontam.controls.out, 
                                           !(sample_control %in% "pos.controls"))
# 8105 taxa in 347 samples
summarize_phyloseq(ps.noncontam.controls.out)

###########################
# prune those not in at least 1 sample
ps.prune <- prune_taxa(taxa_sums(ps.noncontam.controls.out) > 1, ps.noncontam.controls.out) 
# 8056 in 347 samples

# remove samples with < 100 reads
ps.prune <- prune_samples(sample_sums(ps.prune) > 100, ps.prune)
ps.prune
# 8056 taxa in 328 samples = 19 samples lost in pruning


###########################
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps.prune),
               MARGIN = ifelse(taxa_are_rows(ps.prune), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.prune),
                    tax_table(ps.prune))

prev.ASVs<-plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
colnames(prev.ASVs)<- c("Phylum", "mean.prevalence", "sum.prevalence")
write.csv(prev.ASVs, "output/reanalysis/prev.ASVs_local.csv")


rich<-estimate_richness(ps.prune, split = TRUE, measures = NULL)
rownames(rich) <- gsub(x = rownames(rich), pattern = "X", replacement = "")
write.csv(rich, "output/reanalysis/richness.table_local.csv")

########### let's inspect
df.ps.prune <- as.data.frame(sample_data(ps.prune))
df.ps.prune$LibrarySize <- sample_sums(ps.prune) # this is the # of reads
df.ps.prune <- df.ps.prune[order(df.ps.prune$LibrarySize),]
df.ps.prune$Index <- seq(nrow(df.ps.prune))

# library size / number of reads
df.ps.prune$LibrarySize

ps.prune.local<-ps.prune
########### 
saveRDS(ps.prune.local, "output/reanalysis/ps.prune_local.RDS") # save phyloseq
ps.prune.local<-readRDS("output/reanalysis/ps.prune_local.RDS") # bring it back


##### END!!!!
