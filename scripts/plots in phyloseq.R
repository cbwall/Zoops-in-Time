## investigate the final 16S library


remotes::install_github("microbiome/microbiome")
library("microbiome")

PS.fin<-read_csv2phyloseq(
  otu.file = "output/TSCC-output/otu_table.csv",
  taxonomy.file = "output/TSCC-output/tax_table.csv",
  metadata.file = "output/TSCC-output/sam_data.csv",
  sep = ","
)



### plots
# read in formatted metaData
run.metaD<-read.csv("output/TSCC-output/run.metaD.edit.csv")

# format metadata (have to redo this if loading back in)
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Organism", "Time.Point", "Lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

#richness by Lake
richness.plot<-plot_richness(PS.fin, x="Lake", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/richness.plot.location.pdf", height=4, width=10)
dev.off() 

#richness by organism (or sample type i.e, water)
richness.plot<-plot_richness(PS.fin, x="Organism", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/richness.plot.organism.pdf", height=4, width=12)
dev.off() 

#richness by Time.point
richness.plot<-plot_richness(PS.fin, x="Time.Point", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/richness.plot.Time.Point.pdf", height=4, width=10)
dev.off() 


# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(PS.fin))

# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=5000, lty=2)

hist.depth
dev.copy(pdf, "output/hist.depth.pdf", height=4, width=5)
dev.off() 

##
pdf(file="output/rare.raw.pdf", height=4, width=10)
rarecurve(otu_table(PS.fin), step=50, cex=0.5, xlim=c(0,100000), label=FALSE)
#abline(v = 5000, lty = "dotted", col="red", lwd=2)
dev.off() 