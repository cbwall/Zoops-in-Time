## Alpine lake microbiomes: Geography and host identity shape intraseasonal variation of free-living and zooplankton associated microbial communities in alpine lakes.  
<a href="https://doi.org/10.5281/zenodo.15708662"><img src="https://zenodo.org/badge/583512332.svg" alt="DOI"></a>

CB Wall<sup>1,2</sup>, MG Perreault<sup>1</sup>, MY Demmel<sup>1</sup>, EM Diaz<sup>1</sup>, JH Dominguez<sup>1</sup>, JB Shurin<sup>1</sup>. (2025) _Molecular Ecology_.  
<sup>1</sup>University of California, San Diego   
<sup>2</sup>University of Hawai'i at Mānoa  
  
<p align="center">
  <img align="center" src="https://github.com/cbwall/Zoops-in-Time/blob/main/output/Granite downslope.jpg" width="60%">
</p>  
  
  
Microbes contribute to the functio of aquatic ecosystems and zooplankton animals. Many factors affect the taxonomic compositions of free-living (bacterioplankton) and zooplankton-associated microbial communities in lakes, yet how these communities vary seasonally and among lakes remains poorly understood. Here we investigate how free-living bacterial communities and those associated with different crustacean zooplankton hosts change in response to fluctuations in their natural environment across time and space. We repeatedly sampled bacterioplankton, zooplankton communities, and zooplankton microbiomes, and water chemistry parameters of six lakes in the eastern Sierra Nevada mountains of California across the 2022 summer season.  


## File Directory  
The file directory contains folders and scripts (Rmd) to be run in RStudio. The folders house data, output, and raw + polished figures.  
   - `Zoops-in-Time.Rproj` = the R project for the analysis, load this first to allow code to run from correct directory in R Studio
   - `Zoops-in-Time.Rmd` = the scripts and annotated code chunks here will walk through analyses and produce outputs
   - `Zoops-in-Time.html` = the knit output of the Rmd file. Open this in any browser.

   - Folders
     - `data` = contains raw and processed data files
     - `figures` = contains raw exported figures from code in R
     - `output` = contains output subfolders and images 

## Metadata  
Important data files to generate maps, figures, and run models can be found in the *data folder*.  

### Data  
  - `field_nutr_metadata.csv` = field sampling metadata with sampling dates, lake lat-long, lake depth, tow depth, volume of filtered water for baterioplankton (Sterivex), in situ measurments from YSI, and dissolved nutrient analysis.  
  - `data/full.library.env.metaD.csv` = metadata for sequencing run, including: primer sequences, number of zooplankton extracted from, and paired environmental data for site-time sampling.
  - `data/weather.SNARL.csv` = data on rainfall and temperature at the Sierra Nevada Research Laboratory (SNARL) in southern Mammoth Lakes, CA.
  - `data/zoop_community.csv` = sampling metadata with zooplankton community counts (number of organisms), volume of plankton tows, and the relative abundance of plankton.

### Scripts
  - `dada.scripts/` = folder with 4 scripts for dada2 processing in R
  - `cutadapt_16S_git.sh` = cutadapt terminal script for removing primers

### Output  
  - `Granite downslope.jpg` = image from Granite 2 downslope into Yosemite National Park
  - `MetaD.SQ.Env.csv` = dada2 metadata output merged and cleaned used in phyloseq

#### Subfolders 
`Output/` has several subfolders related to microbiome processing.  
  - `dada proc/` = RDS (R data files) from dada and phyloseq.
    - `PS.prune_local` is the final phyloseq object from pipeline after low abundance taxa removed.
    - `PS.fin` the final phyloseq object (unrarified) with merged and re-organized metadata.
    - `PS.500` is rarified phyloseq data to 500 reads per sample.
    - other files include the forward and reverse quality read plots from dada2.
    
  - `final.figures/` = processed/cleaned figures used in publication are in the *output.     
  - `maps/` = images and compiled plots of sample site maps.  
  - `phyloseq output/` = subfolder with files (OTU table, taxonomy, metadata) used to generate phyloseq object (or can be run in other forms).  
  - `phyloseq output_rarified/` = same as above, but 500-read rarified data.  
