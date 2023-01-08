############ 16S part 4
# decontam

# remember to run through the 4.0.2 version of R instead of default


install.packages("ggplot2")
library("ggplot2")

install.packages("decontam")
library("decontam")
packageVersion("decontam")

install.packages("phyloeq")
library("phyloeq")
packageVersion("phyloseq")

install.packages("Biostrings")
library("Biostrings")
packageVersion("Biostrings")

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
print("Tax columns")
colnames(tax_table)


## sample data
# metadata is run.metaD
all(rownames(seqtab.nochim) %in% run.metaD$sampleNames)

rownames(run.metaD) <- run.metaD$sampleNames

ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
              sample_data(run.metaD), 
              tax_table(taxa))


# make a string of DNA names and add to phyloseq
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# save ps object
saveRDS(ps, file="output/ps.rds")

###########################
# Show available ranks in the dataset
rank_names(ps)

table(tax_table(ps)[, "Phylum"], exclude = NULL)
table(tax_table(ps)[, "Kingdom"], exclude = NULL)

# remove NAs in taxonomic table
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Chloroplast"))
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Chloroplast"))

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
saveRDS(ps.prune, file="output/ps.prune")

sample_data(ps.prune)$is.neg <- sample_data(ps.prune)$sample_control == "neg.controls" 
contamdf.prev <- isContaminant(ps.prune, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) # which are contaminants?
head(which(contamdf.prev$contaminant))

###########################
### remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps.prune)
ps.noncontam
saveRDS(ps.noncontam, file="output/ps.noncontam")

rich<-estimate_richness(ps.noncontam, split = TRUE, measures = NULL)
richness.plot<-plot_richness(ps.noncontam, x="year", measures=c("Observed", "Shannon")) + theme_bw()

richness.plot
dev.copy(pdf, "output/richness.plot.pdf", height=4, width=5)
dev.off() 

########### let's inspect
df.noncontam <- as.data.frame(sample_data(ps.noncontam))
df.noncontam$LibrarySize <- sample_sums(ps.noncontam) # this is the # of reads
df.noncontam <- df.noncontam[order(df.noncontam$LibrarySize),]
df.noncontam$Index <- seq(nrow(df.noncontam))
########### 

# library size / number of reads
df.noncontam$LibrarySize

#remove neg controls
ps.noncontam.negout<- subset_samples(ps.noncontam, !(sample_control %in% "neg.controls"))

# 4777 taxa in 99 samples (+ controls still in here -- can use to assess accuracy in mock)

