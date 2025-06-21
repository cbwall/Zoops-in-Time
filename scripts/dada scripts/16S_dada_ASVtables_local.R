############ 16S part 3

# remember to run through the 4.0.2 version of R instead of default
library(dplyr)
packageVersion("dplyr")
library(decontam)
packageVersion("decontam")

start_time <- Sys.time() # track timing

###################################################################################################

#setwd('/projects/ps-shurinlab/users/cbwall/Zoops_MP')

# load sequence table
seqtab.nochim<-readRDS("output/dada proc/seqtab.nochim.rds")

# CREATING ASV FASTA FILE

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# create and write FASTA table with sequences for each ASV
asv_fasta <- c(rbind(asv_headers, asv_seqs))
saveRDS(asv_fasta, "output/dada proc/ASV_fasta.rds")
write(asv_fasta, "output/dada proc/ASVs.fa")


# CREATING COUNTS TABLE
# create counts table by transposing seq table and changing row names
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "output/dada proc/ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)


# CREATING TAXONOMY TABLE
# remove > from headers
asv_headers_tax <- sub('.', '', asv_headers)

taxa <- readRDS("output/dada proc/taxa.rds")

# change row names from sequences to ASV numbers
rownames(taxa) <- asv_headers_tax

write.table(taxa, "output/dada proc/ASV_taxonomy.tsv", sep = "\t", quote = F, col.names = NA)

# final output is a fasta file (and saved R object), counts table, and taxonomy table

