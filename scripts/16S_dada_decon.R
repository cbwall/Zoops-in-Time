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


start_time <- Sys.time() # track timing

setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

## required load in
###################################################################################################
# read in formatted metaData
run.metaD<-read.csv("output/run.metaD.edit.csv")

# Sample Names
sampleNames<-run.metaD$sampleNames

taxa <- readRDS("output/taxa.rds")
seqtab.nochim <- readRDS("output/seqtab.nochim.rds")

###################################################################################################
### end required load in 

# load the three files generated in tables.R
counts <- read.table(file = 'output/ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
tax_table <- read.table(file = 'output/ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
asv_fasta <- readRDS("output/ASV_fasta.rds")

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
# 14094 taxa in 368 samples

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
# 195 Archaea, 12209 Bacteria, 22 Eukaryota 1668 NA

# remove NAs in taxonomic table
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # only Archaea and Bacteria
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Chloroplast"))

# could start with removing ALL non bacterial sequences
# ps<-rm_nonbac(ps)
# by avoiding this step we are keeping in the Archaea

#### summary
# 7577 taxa in 368 samples

###########################
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

###########################
## ID contaminants
# first let's prune those not in at least 1 sample
ps.prune <- prune_taxa(taxa_sums(ps) > 1, ps)

# remove samples with < 100 reads
ps.prune <- prune_samples(sample_sums(ps.prune) > 100, ps.prune)
ps.prune
# 7356 taxa in 341 samples

sample_data(ps.prune)$is.neg <- sample_data(ps.prune)$sample_control == "neg.controls" 
contamdf.prev <- isContaminant(ps.prune, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) # which are contaminants? 7091 NO, 265 YES
head(which(contamdf.prev$contaminant))

###########################
### remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps.prune)
ps.noncontam

#make sure those negative and pos controls are out! 
ps.noncontam.controls.out<- subset_samples(ps.noncontam, 
                                           !(sample_control %in% "neg.controls"))
ps.noncontam.controls.out<- subset_samples(ps.noncontam.controls.out, 
                                           !(sample_control %in% "pos.controls"))

# leaves 322 samples with 7091 taxa


rich<-estimate_richness(ps.noncontam.controls.out, split = TRUE, measures = NULL)
write.csv(rich, "output/richness.table.csv")

########### let's inspect
df.noncontam.noncontrol <- as.data.frame(sample_data(ps.noncontam.controls.out))
df.noncontam.noncontrol$LibrarySize <- sample_sums(ps.noncontam.controls.out) # this is the # of reads
df.noncontam.noncontrol <- df.noncontam.noncontrol[order(df.noncontam.noncontrol$LibrarySize),]
df.noncontam.noncontrol$Index <- seq(nrow(df.noncontam.noncontrol))

# library size / number of reads
df.noncontam$LibrarySize

########### 

# save ps as 4 csvs
# One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv) 
write_phyloseq(ps.noncontam.controls.out, path = "output")



##### END!!!!
