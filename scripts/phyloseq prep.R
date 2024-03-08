## investigate the final 16S library

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

devtools::install_github("benjjneb/dada2", ref="v1.20") # update to most recent dada2
devtools::install_github("zdk123/SpiecEasi")
remotes::install_github("microbiome/microbiome")


# use pacman to load CRAN packages missing
pacman::p_load('knitr', 'microbiome', 'phyloseq', 'tidyr', 'tidyverse', 'knitr', 'magrittr', 'effects', 'devtools',
               'stringi', 'dplyr', "ggplot2", "gridExtra", "dada2", "phyloseq", "vegan", "cowplot", 'doBy', 'ecodist',
               'glue', 'geosphere', 'data.table', 'patchwork', 'car', 'ggcorrplot', 'FactoMineR', 'devtools', 'reshape',
               'lattice',  'plyr', 'magrittr', 'factoextra', 'multcompView', 'decontam', 'factoextra', 'car', 'mia', 
               'ade4', 'fossil', 'picante', 'reshape', 'readr', 'corrr', 'scipplot', 'Hmisc', 
               "decontam","BiocManager", 'ggpubr', 'ggmap', "ggordiplots", "fossil", "SpiecEasi", "igraph", "huge", "MASS")




##########################################################

# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune.LOCAL<-readRDS("output/local/ps.prune_local.RDS")

# idiot check: make sure to remove all mitochondria and chloroplast taxonomic IDs
ps.prune.LOCAL <- subset_taxa(ps.prune.LOCAL, Family!= "Mitochondria" | 
                        is.na(Family) & Class!="Chloroplast" | is.na(Class))

# export metadata and play
metaD<-microbiome::meta(ps.prune.LOCAL)

# adjust caps
metaD<-metaD %>% 
  dplyr::rename(
    sequencing_ID = Sequencing_ID,
    sample_type = Sample_Type,
    organism = Organism,
    time_point= Time.Point,
    lake = Lake,
    sample_ID = Original_ID
  )

# correct NAs to be "water"
metaD$organism<-as.character(metaD$organism) # need to do this to change the level
metaD$organism[is.na(metaD$organism)] = "Water"
metaD$organism<-as.factor(metaD$organism) # back to factor


# format metadata (have to redo this if loading back in)
make.fac<-c("sample_ID", "sequencing_ID", "sample_type", "time_point", "lake", "Plate", "Plate_name", "Well", "sample_control")

metaD[make.fac]<-lapply(metaD[make.fac], factor) # make all these factors

# create phy group (Cladocera vs Copepoda vs water) column
metaD$phy_group <- as.factor(ifelse(metaD$organism == "Calanoid" | 
                            metaD$organism == "Cyclopoid" | 
                            metaD$organism == "Large Calanoid", "copepoda", 
                   ifelse(metaD$organism == "Daphnia" | 
                            metaD$organism == "Ceriodaphnia" | 
                            metaD$organism =="Bosmina" | 
                            metaD$organism == "Holopedium", "cladocera", 
                            "water")))

# create basin column
metaD$basin <- as.factor(ifelse(metaD$lake == "Eastern Brook" | 
                                metaD$lake == "Serene" | 
                                metaD$lake == "Convict", "South", 
                                "North"))


##### Read counts and rarefaction curves
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(ps.prune.LOCAL))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)

# are row names the same? if so, add them in
identical(rownames(sample_sum_df), rownames(metaD))
metaD$read.sum<-sample_sum_df$read.sum

# add in dates, as per metadata
metaD$date<-ifelse(metaD$time_point =="1" | metaD$lake =="Blue", "6/29/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Blue", "7/09/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Blue", "7/21/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Blue", "8/04/2022",
                            ifelse(metaD$time_point =="5" | metaD$lake =="Blue", "8/22/2022",
                    
                     ifelse(metaD$time_point =="1" | metaD$lake =="Convict", "6/25/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Convict", "7/07/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Convict", "7/22/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Convict", "8/05/2022",
                            ifelse(metaD$time_point =="5" | metaD$lake =="Convict", "8/23/2022",
                          
                     ifelse(metaD$time_point =="1" | metaD$lake =="Cooney", "6/29/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Cooney", "7/09/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Cooney", "7/21/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Cooney", "8/04/2022",
                            ifelse(metaD$time_point =="5" | metaD$lake =="Cooney", "8/22/2022", 
                                
                     ifelse(metaD$time_point =="1" | metaD$lake =="Eastern Brook", "6/26/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Eastern Brook", "7/06/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Eastern Brook", "7/22/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Eastern Brook", "8/09/2022",
                            ifelse(metaD$time_point =="5" | metaD$lake =="Eastern Brook", "8/23/2022",
                                   
                    ifelse(metaD$time_point =="1" | metaD$lake =="Serene", "6/26/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Serene", "7/06/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Serene", "7/22/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Serene", "8/09/2022",
                            ifelse(metaD$time_point =="5" | metaD$lake =="Serene", "8/23/2022",
                                  
                    ifelse(metaD$time_point =="1" | metaD$lake =="Virginia", "6/29/2022",
                            ifelse(metaD$time_point =="2" | metaD$lake =="Virginia", "7/09/2022",
                            ifelse(metaD$time_point =="3" | metaD$lake =="Virginia", "7/21/2022",
                            ifelse(metaD$time_point =="4" | metaD$lake =="Virginia", "8/04/2022",
                            "8/22/2022"
                            )))))))))))))))))))))))))))))

# adjust formatting to be Date (from)
metaD$date<-as.Date(mdy(metaD$date))


# create latitude column
metaD$latitude <- as.factor(ifelse(metaD$lake == "Blue", "38.050952", 
                            (ifelse(metaD$lake == "Convict", "37.590094",
                            (ifelse(metaD$lake == "Cooney", "38.04806",
                            (ifelse(metaD$lake == "Eastern Brook", "37.431526", 
                            (ifelse(metaD$lake == "Serene", "37.438392" ,
                            (ifelse(metaD$lake == "Virginia", "38.046975",
                                  "0"))))))))))))

# create longitude column
metaD$longitude <- as.factor(ifelse(metaD$lake == "Blue", "-119.270331", 
                            (ifelse(metaD$lake == "Convict", "-118.857422",
                            (ifelse(metaD$lake == "Cooney", "-119.27702",
                            (ifelse(metaD$lake == "Eastern Brook", "-118.742615", 
                            (ifelse(metaD$lake == "Serene", "-118.744386" ,
                            (ifelse(metaD$lake == "Virginia", "-119.265076",
                                  "0"))))))))))))

metaD<- metaD %>%
  dplyr::select(sequencing_ID, sampleNames, sample_ID, read.sum, time_point, date, lake, latitude, longitude, sample_type, organism, Number.of.Organism, phy_group)

   
# the metadata above is now up to date with info necessary for downstream analysis
# now merge in the environmental data 

################################## ##################################
# read in the environmental data
env.metad<-read.csv("data/full.library.env.metaD.csv")

# make a column for 'Samplenames' and use this as rownames
env.metad$sampleNames<-env.metad$sequencing_ID
env.metad$sampleNames <- gsub("^.{0,3}", "", env.metad$sampleNames)
rownames(env.metad)<- env.metad$sampleNames

# remove the rows not found in the meta-meta data -- these are samples lost in dada2/controls
env.metad<-env.metad[(env.metad$sampleNames %in% metaD$sampleNames),]

env.metad<- env.metad %>%
  dplyr::select(sample_ID, sequencing_ID, sampleNames, sample_type, organism, Number.of.Organism, time_point, lake, 
                elevation_m, temp_C, chla_ug.L, pH, DO_perc, spc, cond, DOC, TDN, TDP)

################################## ##################################
# no merge, subset to make it easier to merge...
env.metad.reduce<-env.metad %>%
  dplyr::select(sampleNames, elevation_m, temp_C, chla_ug.L, pH, DO_perc, spc, cond, DOC, TDN, TDP)


# merge the two = final metadata, add in rownames for phyloseq merge
MetaD.SQ.Env<-merge(metaD, env.metad.reduce, by="sampleNames", all.x=TRUE)
rownames(MetaD.SQ.Env)<-MetaD.SQ.Env$sampleNames

# make new column for above or below treeline
MetaD.SQ.Env$elev.cat <- as.factor(ifelse(MetaD.SQ.Env$elevation_m < 2900, "Below", "Above"))

# are row names the same? if so, let's boogie
identical(rownames(sample_data(ps.prune.LOCAL)), rownames(MetaD.SQ.Env))
write.csv(MetaD.SQ.Env, "output/local/MetaD.SQ.Env.csv")

#updated metadata back into phyloseq
sample_data(ps.prune.LOCAL)<-MetaD.SQ.Env
PS.fin<- ps.prune.LOCAL

########### 
saveRDS(PS.fin, "output/local/PS.fin.RDS") # save phyloseq

########### NMDS

### PCoA
############ PCoA
install.packages("remotes")
remotes::install_github("microbiome/microbiome")

x<-meta(PS.fin) # get metadata from phyloseq


#### ALL
ORD<- ordinate(PS.fin, method='MDS', distance='bray')
NMDS.all<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism)) +
  ggtitle("ALL taxa") +
  theme_classic()   

NMDS.time<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = time_point), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=time_point)) +
  ggtitle("ALL times") +
  theme_classic()  

NMDS.lake<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake)) +
  ggtitle("ALL lakes") +
  theme_classic()  

#####
NMDS.all.tests<-plot_grid(NMDS.all, NMDS.time, NMDS.lake,
                                ncol=3)
NMDS.all.tests
dev.copy(pdf, "output/local/NMDS.all.tests.pdf", height=7, width=16)
dev.off()
#######


### subset data frames
T1<- subset_samples(PS.fin, time_point=="1")
T2<- subset_samples(PS.fin, time_point=="2")
T3<- subset_samples(PS.fin, time_point=="3")
T4<- subset_samples(PS.fin, time_point=="4")
T5<- subset_samples(PS.fin, time_point=="5")

### subset ordinations
ORD.T1 <- ordinate(T1, method='MDS', distance='bray')
ORD.T2 <- ordinate(T2, method='MDS', distance='bray')
ORD.T3 <- ordinate(T3, method='MDS', distance='bray')
ORD.T4 <- ordinate(T4, method='MDS', distance='bray')
ORD.T5 <- ordinate(T5, method='MDS', distance='bray')

### plots
NMDS.ord.T1<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD.T1) +                                                
  geom_point(aes(color = organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism)) +
  ggtitle("Time1") +
  theme_classic()   

NMDS.ord.T1

###
NMDS.ord.T2<-plot_ordination(
  physeq = T2,                                                   
  ordination = ORD.T2) +                                                
  geom_point(aes(color = organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism)) +
  ggtitle("Time2") +
  theme_classic()     

NMDS.ord.T2

####
NMDS.ord.T3<-plot_ordination(
  physeq = T3,                                                   
  ordination = ORD.T3) +                                                
  geom_point(aes(color = organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism))+
  ggtitle("Time3") +
  theme_classic()     

NMDS.ord.T3

###
NMDS.ord.T4<-plot_ordination(
  physeq = T4,                                                   
  ordination = ORD.T4) +                                                
  geom_point(aes(color = organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism))+ 
  ggtitle("Time4") +
  theme_classic()     

NMDS.ord.T4


###
NMDS.ord.T5<-plot_ordination(
  physeq = T5,                                                   
  ordination = ORD.T5) +                                                
  geom_point(aes(color =organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=organism)) + 
  ggtitle("Time5") +
  theme_classic()     

NMDS.ord.T5

#
NMDS.local.tests.org<-plot_grid(NMDS.ord.T1, NMDS.ord.T2, NMDS.ord.T3, NMDS.ord.T4, NMDS.ord.T5,
                                          ncol=3)
NMDS.local.tests.org
dev.copy(pdf, "output/local/NMDS.local.tests.org.pdf", height=9, width=15)
dev.off()



### plots
NMDS.ord.T1<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD.T1) +                                                
  geom_point(aes(color = lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake)) +
  ggtitle("Time1") +
  theme_classic()   

NMDS.ord.T1

###
NMDS.ord.T2<-plot_ordination(
  physeq = T2,                                                   
  ordination = ORD.T2) +                                                
  geom_point(aes(color = lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake)) +
  ggtitle("Time2") +
  theme_classic()     

NMDS.ord.T2

####
NMDS.ord.T3<-plot_ordination(
  physeq = T3,                                                   
  ordination = ORD.T3) +                                                
  geom_point(aes(color = lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake))+
  ggtitle("Time3") +
  theme_classic()     

NMDS.ord.T3

###
NMDS.ord.T4<-plot_ordination(
  physeq = T4,                                                   
  ordination = ORD.T4) +                                                
  geom_point(aes(color = lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake))+ 
  ggtitle("Time4") +
  theme_classic()     

NMDS.ord.T4


###
NMDS.ord.T5<-plot_ordination(
  physeq = T5,                                                   
  ordination = ORD.T5) +                                                
  geom_point(aes(color =lake), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=lake)) + 
  ggtitle("Time5") +
  theme_classic()     

NMDS.ord.T5

#
NMDS.local.tests.org<-plot_grid(NMDS.ord.T1, NMDS.ord.T2, NMDS.ord.T3, NMDS.ord.T4, NMDS.ord.T5,
                                ncol=3)
NMDS.local.tests.org
dev.copy(pdf, "output/local/NMDS.local.tests.LAKE.pdf", height=9, width=15)
dev.off()









################ Now on to plotting

###### some richness plots
#richness by Lake
richness.plot<-plot_richness(PS.fin, x="lake", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/local/richness.plot.location.pdf", height=4, width=10)
dev.off() 

#richness by organism (or sample type i.e, water)
richness.plot<-plot_richness(PS.fin, x="organism", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/local/richness.plot.organism.pdf", height=4, width=12)
dev.off() 

#richness by time_point
richness.plot<-plot_richness(PS.fin, x="time_point", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "output/local/richness.plot.time_point.pdf", height=4, width=10)
dev.off() 
###### 



# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)

hist.depth
dev.copy(pdf, "output/local/hist.depth.pdf", height=4, width=5)
dev.off() 

##
pdf(file="output/local/rare.raw.pdf", height=4, width=10)
rarecurve(otu_table(PS.fin), step=50, cex=0.5, xlim=c(0,100000), label=FALSE)
#abline(v = 5000, lty = "dotted", col="red", lwd=2)
dev.off() 


#### more plots
pdf(file="output/local/read.by.species.pdf", height=4, width=10)
boxplot(MetaD.SQ.Env$read.sum~ MetaD.SQ.Env$organism)
dev.off() 

pdf(file="output/local/read.by.sample.pdf", height=4, width=5)
boxplot(MetaD.SQ.Env$read.sum~ MetaD.SQ.Env$organism)
dev.off() 

pdf(file="output/local/log.reads.sample.pdf", height=4, width=12)
ggplot(MetaD.SQ.Env, aes(x=organism, y=log(read.sum), color=organism)) + geom_boxplot()
dev.off() 

pdf(file="output/local/reads.sample.pdf", height=4, width=7)
ggplot(MetaD.SQ.Env, aes(x=sample_type, y=read.sum, color=organism)) + geom_boxplot()
dev.off() 
##### 

# get out of phyloseq
Shan.rich<-estimate_richness(PS.fin, measures =c("Observed", "Shannon")) # using rarified transformed data
rownames(Shan.rich)<-sapply(str_remove_all(rownames(Shan.rich),"X"),"[") # will remove the X added to the rownames

Shan.rich.df<-merge(data.frame(sample_data(PS.fin)), Shan.rich, by = "row.names") # merge 

### 
shannon.lake<-ggplot(Shan.rich.df, aes(x=lake, y=Shannon)) + 
  geom_boxplot(aes(color=organism)) +
  ggtitle("Shannon richness") + 
  xlab("Lake") +
  theme(axis.title.y = element_blank()) + theme_classic()

shannon.time<-ggplot(Shan.rich.df, aes(x=time_point, y=Shannon)) + 
  geom_boxplot(aes(color=organism)) +
  ggtitle("Shannon richness") + 
  xlab("Time") +
  theme(axis.title.y = element_blank()) + theme_classic()

depth.box<-ggplot(Shan.rich.df, aes(x=organism, y=log(read.sum))) + 
  geom_boxplot(aes(color=organism)) +
  ggtitle("Shannon richness") + 
  xlab("organisms or Water") +
  ylab("log(readsum)") +
  theme(axis.title.y = element_blank()) + theme_classic() +  
  theme(axis.text.x = element_text(size=7))


richness.boxplots<-plot_grid(shannon.lake, shannon.time, depth.box, ncol=3)
richness.boxplots
dev.copy(pdf, "output/local/richness.plots.pdf", height=5, width=20)
dev.off()


