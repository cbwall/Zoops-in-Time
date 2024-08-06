#SRA metadata formatting

run.metaD<- read.csv("data/Zoops_MP_metadata.csv")
run.metaD<-run.metaD %>%
  dplyr::rename("sample_type" = "Sample_Type", 
                "organism" = "Organism", 
                "time_point" = "Time.Point",
                "lake" = "Lake",
                "number_of_organism"= "Number.of.Organism") 

make.fac<-c("Original_ID", "Sequencing_ID", "sample_type", "organism", "time_point", "lake")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors


# ID the - controls
run.metaD$sample_control<-as.factor(run.metaD$sample_type)
run.metaD$sample_control<-gsub(",", "", run.metaD$sample_control) # delete commas in sample names
run.metaD$sample_control[run.metaD$sample_control=='Control PCR_NEG_1' |
                           run.metaD$sample_control=='Control PCR_NEG_2' |
                           run.metaD$sample_control=='Control PCR_NEG_3' |
                           run.metaD$sample_control=='Control PCR_NEG_4' |
                           run.metaD$sample_control=='Control PCR_NEG_5' |
                           run.metaD$sample_control=='Control PCR_NEG_6' |
                           run.metaD$sample_control=='Control PCR_NEG_7' |
                           run.metaD$sample_control=='Control PCR_NEG_8' |
                           run.metaD$sample_control=='Control DNA_CTRL_2' | 
                           run.metaD$sample_control=='Control DNA_CTRL_3' |
                           run.metaD$sample_control=='Control DNA_CTRL_4' |
                           run.metaD$sample_control=='Control DNA_CTRL_5' |
                           run.metaD$sample_control=='Control DNA_CTRL_6' | 
                           run.metaD$sample_control=='Control DNA_CTRL_7' |
                           run.metaD$sample_control=='Control WATER_DNA_CTRL_1' |
                           run.metaD$sample_control=='Control WATER_DNA_CTRL_2' |
                           run.metaD$sample_control=='DNA_CTRL_1'] <- "neg.controls"

run.metaD$sample_control[run.metaD$sample_control=='Control PCR_POS_1' |
                           run.metaD$sample_control=='Control PCR_POS_2' |
                           run.metaD$sample_control=='Control PCR_POS_3' |
                           run.metaD$sample_control=='Control PCR_POS_4'] <- "pos.controls"


run.metaD$sample_control<-recode_factor(run.metaD$sample_control,
                                        "Zooplankton" ="zooplankton",
                                        "Water" = "water") 

# load in the cut-adapt samples in the "trimmed" folder
# run this in terminal 'gunzip data/trimmed_sequences/*'
# now unzipped... proceed

# path to folder containing demultiplexed library sequencing files
# make sure unzipped
miseq_path<-"data/trimmed" 
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

Sequencing_ID<-paste("MP",sampleNames) # compile
Sequencing_ID<-gsub(" ", "_", Sequencing_ID) # remove space and add an underscore

# remove the "Undetermined" sample name 'S0_R1'
sampleNames<-sampleNames[1:368]

SN.FQ<-as.data.frame(cbind(Sequencing_ID,sampleNames))

# if need to merge by a column, say if sequences not all in a single run or separated for some reason...
run.metaD.merge <- merge(run.metaD, SN.FQ, by="Sequencing_ID")

####### other metadata additions
metaD<-run.metaD.merge
metaD$date_simple<-recode_factor(metaD$time_point, 
                                 "1" = "25-Jun-2022",
                                 "2" = "6-Jul-2022",
                                 "3" = "21-Jul-2022",
                                 "4" = "4-Aug-2022",
                                 "5" = "22-Aug-2022")

metaD$organism<-as.character(metaD$organism)
metaD$organism<-ifelse(metaD$sample_type=="Water", "water", metaD$organism)
metaD$organism[is.na(metaD$organism)] = "controls"

metaD$organism<-as.factor(metaD$organism)
metaD$organism<-factor(metaD$organism, levels =c("Bosmina", "Ceriodaphnia", "Daphnia", "Holopedium", 
                                                 "Cyclopoid", "Calanoid", "Large Calanoid",
                                                 "water", "controls")) # back to factor

metaD$organism<-recode_factor(metaD$organism, "Large Calanoid" = "Calanoid")


# create phy group (Cladocera vs Copepoda vs water) column
metaD$phy_group <- as.factor(ifelse(metaD$organism == "Calanoid" | 
                                      metaD$organism == "Cyclopoid", "copepoda", 
                                    ifelse(metaD$organism == "Daphnia" | 
                                             metaD$organism == "Ceriodaphnia" | 
                                             metaD$organism =="Bosmina" | 
                                             metaD$organism == "Holopedium", "cladocera",
                                           ifelse(metaD$organism == "water", "water", "controls"))))

# add in dates, as per metadata
metaD$collection_date<-ifelse(metaD$time_point =="1" | metaD$lake =="Blue", "6/29/2022",
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


# create basin column
metaD$basin <- as.factor(ifelse(metaD$lake == "Eastern Brook" | 
                                  metaD$lake == "Serene" | 
                                  metaD$lake == "Convict", "South", 
                                "North"))

# create latitude column
metaD$latitude <- as.numeric(ifelse(metaD$lake == "Blue", "38.050952", 
                                    (ifelse(metaD$lake == "Convict", "37.590094",
                                            (ifelse(metaD$lake == "Cooney", "38.04806",
                                                    (ifelse(metaD$lake == "Eastern Brook", "37.431526", 
                                                            (ifelse(metaD$lake == "Serene", "37.438392" ,
                                                                    (ifelse(metaD$lake == "Virginia", "38.046975",
                                                                            "0"))))))))))))

# create longitude column
metaD$longitude <- as.numeric(ifelse(metaD$lake == "Blue", "-119.270331", 
                                     (ifelse(metaD$lake == "Convict", "-118.857422",
                                             (ifelse(metaD$lake == "Cooney", "-119.27702",
                                                     (ifelse(metaD$lake == "Eastern Brook", "-118.742615", 
                                                             (ifelse(metaD$lake == "Serene", "-118.744386" ,
                                                                     (ifelse(metaD$lake == "Virginia", "-119.265076",
                                                                             "0"))))))))))))

# change to N and W for lat long
metaD$latitude <- paste0(metaD$latitude, " N")
metaD$latitude<- ifelse(metaD$latitude=="NA N", NA, metaD$latitude)

# make +
metaD$longitude<-abs(metaD$longitude)
metaD$longitude <- paste0(metaD$longitude, " W")
metaD$longitude<- ifelse(metaD$longitude=="NA W", NA, metaD$longitude)


metaD$lat_long<-interaction(metaD$latitude, metaD$longitude, sep=" ")
metaD$sequences<-"bacteria"
metaD$depth<-"NA"
metaD$env_broad_scale<-"aquatic microbiome"
metaD$env_local_scale<- "alpine lake microbiome"
metaD$env_medium<- ifelse(metaD$sample_control=="zooplankton", "zooplankton microbiome", 
                               ifelse(metaD$sample_control=="water", "bacterioplankton", "controls"))
metaD$geo_loc_name<-"USA:California,Eastern_Sierra_mountains"

# rename factors
metaD <- metaD %>% 
  dplyr::rename("taxa_water_control" = "organism",
                "sample_names_short" = "sampleNames")

  
metaD.SRA<-metaD %>%
  dplyr::select(Sequencing_ID, sequences, collection_date, depth, env_broad_scale, 
                env_local_scale, env_medium, geo_loc_name, lat_long, time_point, date_simple, collection_date, lake, 
                basin, sample_control, phy_group, taxa_water_control, number_of_organism, sample_names_short)

write.csv(metaD.SRA, "output/reanalysis/run.metaD.SRA.csv")


#############################
##########################################################

# samplenames
# path to folder containing demultiplexed library sequencing files
# make sure unzipped
fastq_path<-"data/fastq" 

#full file names
fq.names<-list.files(fastq_path, ".gz")

#sequence names, F/R
fq.names.noext<-sub('\\_001.fastq.gz$', '', fq.names) 

### Sort forward/reverse reads
# these are forward files
f.fastq <- sort(list.files(fastq_path, pattern="_R1_001.fastq.gz"))

# these are reverse files
r.fastq <- sort(list.files(fastq_path, pattern="_R2_001.fastq.gz"))

##
fq.Names.p1 <- sapply(strsplit(f.fastq, "_"), `[`, 1) # extract sample names
fq.Names.p2 <- sapply(strsplit(f.fastq, "_"), `[`, 2) # extract sample names
fq.Names.p3 <- sapply(strsplit(f.fastq, "_"), `[`, 3) # extract the run # sample
fq.Names.p4 <- sapply(strsplit(f.fastq, "_"), `[`, 4) # extract the run # sample

fq.Names.short<-paste(fq.Names.p1, fq.Names.p2, fq.Names.p3) # compile
fq.Names.short<-gsub(" ", "_", fq.Names.short) # remove space and add an underscore

fq.Names.long<-paste(fq.Names.p1, fq.Names.p2, fq.Names.p3, fq.Names.p4) # compile
fq.Names.long<-gsub(" ", "_", fq.Names.long) # remove space and add an underscore

FQ.SRA<-as.data.frame(cbind(fq.Names.short, fq.Names.long, f.fastq, r.fastq))
write.csv(FQ.SRA, "output/reanalysis/FQ.SRA.csv")
