## investigate the final 16S library


remotes::install_github("microbiome/microbiome")
library("microbiome")
library("vegan")
library("tidyr")

# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune<-readRDS("output/ps.prune.RDS")

# read in metadata and re-format so it has the str necessary
M<-read.csv("output/sam_data.csv")

# set row names
row.names(M)<-M$X 

# remove junk column
M<-M[-1]

# format metadata (have to redo this if loading back in)
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Time.Point", "Lake", "Plate", "Plate_name", "Well", "sample_control")

M[make.fac]<-lapply(M[make.fac], factor) # make all these factors

M$Organism<-M$Organism %>% replace_na("Water")
M$Organism<- as.factor(M$Organism)


#### assemble the phyloseq object from raw files
PS.fin<-
  read_phyloseq(
    otu.file = "output/otu_table.csv",
    taxonomy.file = "output/tax_table.csv",
    metadata.file = "output/sam_data.csv",
    sep = ","
  )

#replace sample data with reformatted metadata
sample_data(PS.fin)<-M

# check that levels exist
levels(get_variable(PS.fin, "Lake"))
levels(get_variable(PS.fin, "Organism"))

################ Now on to plotting

###### some richness plots
#richness by Lake
richness.plot<-plot_richness(PS.fin, x="Lake", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.location.pdf", height=4, width=10)
dev.off() 

#richness by organism (or sample type i.e, water)
richness.plot<-plot_richness(PS.fin, x="Organism", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.organism.pdf", height=4, width=12)
dev.off() 

#richness by Time.point
richness.plot<-plot_richness(PS.fin, x="Time.Point", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.Time.Point.pdf", height=4, width=10)
dev.off() 
###### 

##### Read counts and rarefaction curves
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(PS.fin))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)

# merge in the reads
run.metaD<-merge(M, sample_sum_df, by="sampleNames", all.y=TRUE)
write.csv(run.metaD, "output/run.metaD.final.csv")


# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)

hist.depth
dev.copy(pdf, "figures/hist.depth.pdf", height=4, width=5)
dev.off() 

##
pdf(file="figures/rare.raw.pdf", height=4, width=10)
rarecurve(otu_table(PS.fin), step=50, cex=0.5, xlim=c(0,100000), label=FALSE)
#abline(v = 5000, lty = "dotted", col="red", lwd=2)
dev.off() 


#### more plots
pdf(file="figures/read.by.species.pdf", height=4, width=10)
boxplot(run.metaD$read.sum~run.metaD$Organism)
dev.off() 

pdf(file="figures/read.by.sample.pdf", height=4, width=5)
boxplot(run.metaD$read.sum~run.metaD$sample_control)
dev.off() 

pdf(file="figures/log.reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=log(read.sum), color=Organism)) + geom_boxplot()
dev.off() 

pdf(file="figures/reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=read.sum, color=Organism)) + geom_boxplot()
dev.off() 
##### 





### PCoA
############ PCoA
sample_variables(PS.fin)

# make colors for sites
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


###
T1<- subset_samples(PS.fin, Time.Point=="1")
ORD.T1 <- ordinate(T1, method='MDS', distance='bray')

NMDS.ord.T1<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD.T1) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.T1
dev.copy(pdf, "figures/NMDS.ord.T1.pdf", height=6, width=7)
dev.off() 

###
T2<- subset_samples(PS.fin, Time.Point=="2")
ORD.T2 <- ordinate(T2, method='MDS', distance='bray')

NMDS.ord.T2<-plot_ordination(
  physeq = T2,                                                   
  ordination = ORD.T2) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time2") +
  theme_classic()     

NMDS.ord.T2
dev.copy(pdf, "figures/NMDS.ord.T2.pdf", height=6, width=7)
dev.off() 

####
T3<- subset_samples(PS.fin, Time.Point=="3")
ORD.T3 <- ordinate(T3, method='MDS', distance='bray')

NMDS.ord.T3<-plot_ordination(
  physeq = T3,                                                   
  ordination = ORD.T3) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_shape_manual(values=c(19,17,15,3,7,8,2,0)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time3") +
  theme_classic()     

NMDS.ord.T3
dev.copy(pdf, "figures/NMDS.ord.T3.pdf", height=6, width=7)
dev.off() 


###
T4<- subset_samples(PS.fin, Time.Point=="4")
ORD.T4 <- ordinate(T4, method='MDS', distance='bray')

NMDS.ord.T4<-plot_ordination(
  physeq = T4,                                                   
  ordination = ORD.T4) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time4") +
  theme_classic()     

NMDS.ord.T4
dev.copy(pdf, "figures/NMDS.ord.T4.pdf", height=6, width=7)
dev.off() 


###
T5<- subset_samples(PS.fin, Time.Point=="5")
ORD.T5 <- ordinate(T5, method='MDS', distance='bray')

NMDS.ord.T5<-plot_ordination(
  physeq = T5,                                                   
  ordination = ORD.T5) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time5") +
  theme_classic()     

NMDS.ord.T5
dev.copy(pdf, "figures/NMDS.ord.T5.pdf", height=6, width=7)
dev.off() 

