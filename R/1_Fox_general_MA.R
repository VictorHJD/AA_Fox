## Please un comment the first time you run this and re-install packages

## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

## library(MultiAmplicon)

## using the dev version!
devtools::load_all("../MultiAmplicon")

library(ggplot2)
library(dada2)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(tidyr)
library(dplyr)
library(stargazer)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE

doMultiAmp <- FALSE

doTax <- FALSE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run foxes#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/2018_22_fox_all/" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[205]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[100]])

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/filtered_complete"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(250,250), minLen=c(250,250), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file 

#Primers used in the arrays 
ptable <- read.csv(file = "input_data/primer_file_foxes.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)
MAF <- MultiAmplicon(primer, files)

##Multi amplicon pipeline
if(doMultiAmp){
    filedir <- "/SAN/Metabarcoding/AA_Fox/stratified_files/"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MAF <- sortAmplicons(MAF, filedir=filedir)
  
    errF <-  learnErrors(unlist(getStratifiedFilesF(MAF)), 
                         verbose=0, multithread = 48)
    errR <- learnErrors(unlist(getStratifiedFilesR(MAF)), 
                        verbose=0, multithread = 48)
  
    MAF <- dadaMulti(MAF, Ferr=errF, Rerr=errR,  pool=FALSE,
                     verbose=0, mc.cores=48)

    MAF <- mergeMulti(MAF, mc.cores=48) 
  
    propMerged <- MultiAmplicon::calcPropMerged(MAF)
  
    summary(propMerged)
    table(propMerged<0.8)
  
    MAF <- mergeMulti(MAF, justConcatenate=propMerged<0.8, mc.cores=48) 
  
    MAF <- makeSequenceTableMulti(MAF, mc.cores=48)
    

    ## fill it, bind it, coerce it to integer
    STF <- getSequenceTable(MAF, dropEmpty=FALSE)
    STFU <- do.call(cbind, STF)
    mode(STFU) <- "integer"
    
    STFU <- STFU[, !duplicated(colnames(STFU))]
    ## and very very harsh chimera removal
    isCruelBimera <- dada2::isBimeraDenovoTable(STFU, multithread=TRUE,
                                                minSampleFraction=0.5,
                                                allowOneOff=TRUE, maxShift = 32,
                                                ignoreNNegatives=4)

    ## the same pooled over all samples
    isPooledBimera <- dada2::isBimeraDenovo(STFU, multithread=TRUE, 
                                            allowOneOff=TRUE, maxShift = 32)

    ## saveRDS(isCruelBimera, "/SAN/Metabarcoding/AA_Fox/cruelBimera.Rds")
    ## saveRDS(isPooledBimera, "/SAN/Metabarcoding/AA_Fox/pooledBimera.Rds")

    table(isCruelBimera ,  isPooledBimera)
    
    superCruel <- isCruelBimera | isPooledBimera

    sum(STFU[, !superCruel])
    sum(STFU[, !isCruelBimera])
    sum(STFU[, !isPooledBimera])
    
    NoBimeras <- colnames(STFU[, !superCruel])
  
    foo <- lapply(getSequenceTable(MAF), function (x) {
        noBim <- intersect(NoBimeras, colnames(x))
        x[, noBim]
    })
    MAF@sequenceTableNoChime <- foo
    
    saveRDS(MAF, "/SAN/Metabarcoding/AA_Fox/MAF_complete.RDS")
} else{
    MAF <- readRDS("/SAN/Metabarcoding/AA_Fox/MAF_complete.RDS") 
}


###  TO FIX IN PACKAGE
## tracking <- getPipelineSummaryX(MAF)
## plotPipelineSummary(tracking)

png("figures/suppl/AmpSampleHeatmapRAW.png", width=24, height=8, units = 'in', res = 300)
sumPheatmap <- plotAmpliconNumbers(MAF) ### 
dev.off()

## everything clustering with Negative controls should be excluded!!
SampleClusters <- cutree(sumPheatmap$tree_col, 2)

## but we exclude later
exclude.samples <- names(SampleClusters)[SampleClusters==2]
exclude.samples <- colnames(MAF)%in%exclude.samples


###New taxonomic assignment 
if (doTax){ ## simply save the blast files, that's even faster than
    unlink("/SAN/Metabarcoding/AA_Fox/in.fasta")
    unlink("/SAN/Metabarcoding/AA_Fox/out.blt")
}


## setting doTax to FALSE and re-loading the object
MAF2 <- blastTaxAnnot(MAF,  
                      negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                      db = "/SAN/db/blastdb/nt/nt",
                      infasta = "/SAN/Metabarcoding/AA_Fox/in.fasta",
                      outblast = "/SAN/Metabarcoding/AA_Fox/out.blt",
                      num_threads = 96)
##to phyloseq

### PACKAGE PROBLEM WITH NULL SLOTs
## exclude.primers <- unlist(lapply(getSequenceTableNoChime(MAF2),
##                          function (x) is.null(dim(x))))

## but also problem with empty tables
exclude.primers <- unlist(lapply(getSequenceTableNoChime(MAF2),
                         function (x) {
                             dix <- dim(x)
                             is.null(dix)|!all(dix>0)
                         }
                         ))

### NAME selection doesn't work PACKAGE!!!
## We exclude like this
MAFinal <- MAF2[which(!exclude.primers), which(!exclude.samples)]

png("figures/suppl/AmpSampleHeatmap.png", width=24, height=8, units = 'in', res = 300)
plotAmpliconNumbers(MAFinal)
dev.off()

PS.l <- toPhyloseq(MAFinal,
                   samples=colnames(MAFinal),
                   multi2Single=FALSE)

### how many of these have Nematodes and Platyhelminthes
lapply(PS.l, function(x) {
    any(tax_table(x)[,"phylum"]%in%c("Nematoda", "Platyhelminthes"))
}) %>% unlist() %>% table()


PS <- toPhyloseq(MAFinal,
                 samples=colnames(MAFinal),
                 multi2Single=TRUE)

##Add real sample data
sample.data <- readRDS("intermediate_data/Fox_data_envir.RDS")

##Date for a season category
sample.data$date <- as.Date(sample.data$date_found, "%d.%m.%Y")

## this regular seasons would mean only 5 summer samples... 
table(ifelse(month(sample.data$date)>11, "winter",
      ifelse(month(sample.data$date)>8, "autumn",
      ifelse(month(sample.data$date)>5, "summer",
      ifelse(month(sample.data$date)>2, "spring", "winter")))))

## autumn spring summer winter 
##     69     41      5     99 

## let's  spread them to spring and autumn
table(ifelse(month(sample.data$date)>11, "winter",
      ifelse(month(sample.data$date)>6, "S_autumn",
      ifelse(month(sample.data$date)>2, "spring", "winter"))))


sample.data$season <- ifelse(month(sample.data$date)>11, "winter",
                      ifelse(month(sample.data$date)>6, "S_autumn",
                      ifelse(month(sample.data$date)>2, "spring", "winter")))

sample.data$year <- year(sample.data$date)


##### IMPUTATION of NA missing values #################

### iputing the season NAs by consecutive sampling reasoning, there
### are some cases in which order put those in between seasons. We
### simply ignor and use the season of the sample point later (as they
### came sometimes in batches).
sample.data <- sample.data[order(as.numeric(sample.data$IZW_ID)), ]

table(is.na(sample.data$season))
## 16 missing dates and hence missing seasons

sample.data$SY_imputed <- ifelse(is.na(sample.data$season), TRUE, FALSE)

sample.data %>% fill(season, year, .direction = "updown") ->
    sample.data

### imputing missing weight (after making it numeric
sample.data$weight_kg <- as.numeric((gsub(" *", "", sample.data$weight_kg)))

table(is.na(sample.data$weight_kg))
## 8 missing values for weight

## using the mean of the age, sex and condition
sample.data %>%
    group_by(age, condition, sex) %>%
    mutate(weight_imputed = ifelse(is.na(weight_kg), TRUE, FALSE)) %>%
    mutate(weight_kg = case_when(is.na(weight_kg) ~ mean(weight_kg, na.rm=TRUE),
                                 TRUE ~ as.numeric(weight_kg))) ->
    sample.data

## imputing missing condition
table(is.na(sample.data$condition))

## we use "excellent" if above group mean, "autolytic" if below
sample.data %>%
    group_by(age, sex) %>%
    mutate(condition_imputed = ifelse(is.na(condition), TRUE, FALSE)) %>%
    mutate(condition = case_when(is.na(condition) ~
                                     ifelse(weight_kg > mean(weight_kg, na.rm=TRUE),
                                            "excellent", "autolytic"),
                                 TRUE ~ condition)) ->
    sample.data


sample.data$IZW_ID <- as.vector(sample.data$IZW_ID)

### We have to fix the IZW sample names!
setdiff(sample_names(PS), sample.data$IZW_ID)
## they have these additional "b"s

## just add "b"s to all those
sample.data$IZW_ID[!sample.data$IZW_ID%in%sample_names(PS)] <-
    paste0(sample.data$IZW_ID[!sample.data$IZW_ID%in%sample_names(PS)], "b")

### Now in the other direction
setdiff(sample.data$IZW_ID, sample_names(PS))
## they have these additional "b"s

rownames(sample.data) <- sample.data$IZW_ID

## adding "Fox" for better names
rownames(sample.data) <- paste("Fox", rownames(sample.data))
rownames(PS@sam_data) <- paste("Fox", rownames(PS@sam_data))
sample_names(PS) <- paste("Fox", sample_names(PS))

## align and cbind to get the combinded sample data
PS@sam_data <- sample_data(cbind(PS@sam_data, sample.data[rownames(sample_data(PS)), ]))

## make weight numeric
PS@sam_data[, "weight_kg"] <- as.numeric(unlist(PS@sam_data[, "weight_kg"] ))


################ Analysing DNA quality and quantity ###########
### ############ and addint it to the sample data ########

### we do this without the "b" post-fixes in sample names
sample.data$ID <- paste("Fox", gsub(" |b", "", sample.data$IZW_ID))

## The sampples producing good sequencing data ###################
SD <- sample_data(PS)
SD$ID <- gsub("b", "", rownames(SD))


## All obtained sampels #########################################
Meta <- read.csv("input_data/AllSamples.csv")
Meta <- Meta[!is.na(Meta$IZW.ID), ]
Meta$ID <- paste("Fox", gsub(" ", "", Meta$IZW.ID))

length(unique(Meta$ID))
## 619 foxes sampled                               <-- RESULT?!
dupes <- Meta$ID[duplicated(Meta$ID)]
Meta[Meta$ID%in%dupes, ]
## strange errors in raw data for three samples

### but lukiely these samples were not used furhter
table(SD$ID%in%dupes)
## simply ignore them
Meta <- Meta[!duplicated(Meta$ID), ]

## And their DNA measurements
DNA <- read.csv("input_data/nanodrop_1-4_zusammen.csv")
DNA$ID <- paste("Fox", gsub(" |b|\\.2", "", DNA$Sample.ID))

length(DNA$ID)
length(unique(DNA$ID)) ## 474 foxes processed

## for those DNA was extracted twice we use the sample with the higher
## DNA quantity
DNA %>% group_by(ID) %>%
    filter(ng.ul == max(ng.ul)) -> DNA

length(DNA$ID) ## unique now!!
length(unique(DNA$ID)) ## 474 foxes processed

## merge the Metadata and the DNA data
### 473 match                               <<- RESULT, extracted DNA!!
DNA <- merge(DNA, Meta, by="ID")

length(DNA$ID) ## unique still!!
length(unique(DNA$ID)) ## 474 foxes processed

DNA$sequenced <- DNA$ID%in%sample.data$ID
table(DNA$sequenced)
## 226                                       <<- RESULT, sequenced!

DNA$analysed <- DNA$ID%in%SD$ID
table(DNA$analysed) ## 152, so we're missing DNA qual for 3 foxes

DNA$age..j.a.[DNA$age..j.a.%in%"juvenil"] <- "juvenile"

colnames(DNA) <- gsub("\\..*", "", colnames(DNA))

DNA$weight_kg <- as.numeric(DNA$weight)

colnames(DNA)[colnames(DNA)%in%"ng"] <- "DNAng.ul"

colnames(DNA)[colnames(DNA)%in%"X260"] <- c("DNA260.280", "DNA260.230")


#### confirming: sequenced were the better samples (in terms of quantity)
tapply(DNA$DNAng.ul, DNA$sequenced, median)
tapply(DNA$DNAng.ul, DNA$sequenced, mean)

### 260.230 ratio was important for selection for sequencing 
tapply(DNA$DNA260.230, DNA$sequenced, median)
tapply(DNA$DNA260.230, DNA$sequenced, mean)

### 260.280 ratio was NOT important for selection for sequencing 
tapply(DNA$DNA260.280, DNA$sequenced, median)
tapply(DNA$DNA260.280, DNA$sequenced, mean)



## lower DNA quantity sampels tend to not be analysed futher... they
## failed in sequencing and were screened
tapply(DNA[DNA$sequenced, "DNAng.ul"], DNA[DNA$sequenced, "analysed"], median)
tapply(DNA[DNA$sequenced, "DNAng.ul"], DNA[DNA$sequenced, "analysed"], mean)

## again the same for DNA quality in 260.230 tend to not be analysed
## futher... they failed in sequencing and were screened
tapply(DNA[DNA$sequenced, "DNA260.230"], DNA[DNA$sequenced, "analysed"], median)
tapply(DNA[DNA$sequenced, "DNA260.230"], DNA[DNA$sequenced, "analysed"], mean)

## and not so for 260.280
tapply(DNA[DNA$sequenced, "DNA260.280"], DNA[DNA$sequenced, "analysed"], median)
tapply(DNA[DNA$sequenced, "DNA260.280"], DNA[DNA$sequenced, "analysed"], mean)


############## OVERALL EXTRACTION #############################################

modQuantExtr <- lm(DNAng.ul ~ weight_kg + condition + sex + age + area,
                   data=DNA)
summary(modQuantExtr)
### lower quantity of DNA for juvenile and heavier foxes, much lower
### in Brandenburg!!!

modQual1Extr <- lm(DNA260.280 ~ weight_kg + condition + sex + age + area,
                   data=DNA)
summary(modQual1Extr)
### different (lower) quality in Brandenburg

modQual2Extr <- lm(DNA260.230 ~ weight_kg + condition + sex + age + area,
    data=DNA)
summary(modQual2Extr)
### different (lower) quality in Brandenburg

############## WAS BIASED #####################################################


############## SEQUENCED SAMPLES  #############################################
modQuantSEQ <- lm(DNAng.ul ~ weight_kg + condition + sex + age + area,
                  data=subset(DNA, DNA$sequenced))
summary(modQuantSEQ)
### lower quantity of DNA for juvenile and heavier foxes, much lower
### in Brandenburg!!!

modQual1SEQ <- lm(DNA260.280 ~ weight_kg + condition + sex + age + area,
                   data=subset(DNA, DNA$sequenced))
summary(modQual1Extr)
### different (lower) quality in Brandenburg

modQual2SEQ <- lm(DNA260.230 ~ weight_kg + condition + sex + age + area,
                  data=subset(DNA, DNA$sequenced))
summary(modQual2SEQ)
### different (lower) quality in Brandenburg

############## STILL BIASED #####################################################

############## ANALYSED SAMPLES  #############################################
modQuantANA <- lm(DNAng.ul ~ weight_kg + condition + sex + age + area,
                  data=subset(DNA, DNA$analysed))
summary(modQuantANA)
### lower quantity of DNA for juvenile and heavier foxes and lower
### in Brandenburg... significance gone... but still worriesome

modQual1ANA <- lm(DNA260.280 ~ weight_kg + condition + sex + age + area,
                   data=subset(DNA, DNA$analysed))
summary(modQual1ANA)
### different (lower) quality in Brandenburg

modQual2ANA <- lm(DNA260.230 ~ weight_kg + condition + sex + age + area,
                  data=subset(DNA, DNA$analysed))
summary(modQual2ANA)
### different (lower) quality in Brandenburg

stargazer(modQuantExtr, modQuantSEQ, modQuantANA, type="html", out="tables/suppl/quant.html")
stargazer(modQual1Extr, modQual1SEQ, modQual1ANA, type="html", out="tables/suppl/qual1.html")
stargazer(modQual2Extr, modQual2SEQ, modQual2ANA, type="html", out="tables/suppl/qual2.html")

############## BIAS is REMOVED only for DNA quantity in processing and
############## potentially only partially !!!

PS@sam_data$IDb <- rownames(PS@sam_data)
PS@sam_data$ID <- gsub("b", "", rownames(PS@sam_data))

Sdat <- as.data.frame(as(sample_data(PS), "matrix"))

DNAData <- subset(DNA, DNA$analysed)[, c("ID", "DNAng.ul", "DNA260.230", "DNA260.280")]

Sdat <- merge(Sdat, DNAData, by = "ID", all=TRUE)
rownames(Sdat) <- Sdat$IDb


## they are still alinged
all(rownames(Sdat) == rownames(PS@sam_data))
PS@sam_data <- sample_data(Sdat)

###
## Store all this in the central object of the analysis/repository
## (for reproducibilty)
saveRDS(PS, file="intermediate_data/PhyloSeqCombi.Rds")

###For primer analysis (Victor), still stored on our server 
## saveRDS(PS.l, file="/SAN/Metabarcoding/AA_Fox/PhyloSeqList.Rds") 

## For Caro and the paper, previously on the server, 
## saveRDS(PS, file="/SAN/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds")


### and we make the SYNONYMES Capillaria aerophila and Eucoleus
### aerophilus the same thing

table(tax_table(PS)[tax_table(PS)[, "genus"]%in%
                    "Eucoleus", "species"])  

table(tax_table(PS)[tax_table(PS)[, "genus"]%in%
                    "Capillaria", "species"])
    
tax_table(PS)[tax_table(PS)[, "genus"]%in%
                    "Capillaria", "genus"] <- "Eucoleus"

tax_table(PS)[tax_table(PS)[, "genus"]%in%
                    "Eucoleus", "species"] <- "Eucoleus aerophilus"

table(tax_table(PS)[tax_table(PS)[, "genus"]%in%
                    "Eucoleus", "family"])

### Adding Categories to the taxonomy information!
### We have to addd this to genus agglommerted data!
## collapse to genus level
PSG <- phyloseq::tax_glom(PS, "genus")

## We take the helminth traits information from manually curated our
## input data
traits <- read.csv("input_data/helminth_traits.csv")

traits %>%
    tibble::column_to_rownames("t.genus")  -> traits 

BADtaxa <- rownames(traits)[!traits$BlastEvaluation%in%"Okay"]

NOPara <- rownames(traits)[traits$fox.parasite%in%"No"]

## other non-fox parasites
OtherPara <- rownames(traits)[traits$fox.parasite%in%"No" &
                              !traits$lifecycle%in%"free.living"]

## remove bad taxa
## only 3784 reads for bad taxa when excluding only the bad blast annotation
sum(otu_table(subset_taxa(PSG, genus%in%BADtaxa)))

## only 7842 reads for non-parasitic taxa ... 9939 with the new
## samples previously lost in tables...
sum(otu_table(subset_taxa(PSG, genus%in%NOPara)))

## ## results reporting
PS
PSG
subset_taxa(PS, phylum %in% c("Nematoda", "Platyhelminthes"))
subset_taxa(PSG, phylum %in% c("Nematoda", "Platyhelminthes"))

## Store the ASV naming in the taxtable itself (instead of only in the
## rownames)
tax_table(PSG) <- cbind(tax_table(PSG), ASVn=rownames(tax_table(PSG)))

foo <- merge(tax_table(PSG), traits,
             by.x="genus", by.y=0, all.x=TRUE)

## let's use the genus name as rowname
rownames(foo) <- foo$ASVn

## add the "category" information to the taxonomy

foo$category <-
    ifelse(foo$phylum %in% c("Nematoda", "Platyhelminthes") &
           foo$fox.parasite %in%"Yes", "Helminth", ### <- this is a helminth!
    ifelse((foo$phylum %in% c("Nematoda", "Platyhelminthes") &
            foo$fox.parasite %in%"No") |
           (foo$phylum %in% c("Annelida", "Arthropoda", "Chordata",
                              "Mollusca", "Streptophyta")&
            !foo$genus%in%c("Vulpes", "Homo", "Globicephala",
                            "Hylobates", "Procyon", "Canis")), 
           "Diet", ## <- this is Diet
    ifelse(foo$phylum%in%c("Actinobacteria", "Bacteroidetes",
                           "Deferribacteres", "Firmicutes", "Fusobacteria",
                           "Proteobacteria", "Spirochaetes", "Tenericutes"),
           "Microbiome", ## <- this is the bacterial microbiome
    ifelse(foo$order%in%"Eucoccidiorida",
           "ApicoParasites", ## <- these are apicomplexan parasites
    ifelse(foo$"phylum"%in% c("Ascomycota", "Basidiomycota", "Blastocladiomycota",
                              "Chytridiomycota", "Cryptomycota", "Mucoromycota",
                              "Zoopagomycota"),
           "FungalMicrobiome",  ## <- this is the fungal microbiome
           "noClue" ## <- undecided about everything else!!
    )))))


tax_table(PSG) <- as.matrix(foo)[rownames(tax_table(PSG)), ]

## check that taxa were not messed up and we are still removing the
## bad ones
sum(otu_table(subset_taxa(PSG, genus%in%BADtaxa)))
sum(otu_table(subset_taxa(PSG, genus%in%NOPara)))

## Finally drop the "bad stuff" from our phyloseq object
### Exclude the bad taxa
PSG <- subset_taxa(PSG, !genus%in%BADtaxa)

## the categories of gut content
table(tax_table(PSG)[, "category"])

### Short overview of what has been imputed in the final dataset
## season/year
table(sample_data(PSG)$SY_imputed)

## weight
table(sample_data(PSG)$weight_imputed)

## condition
table(sample_data(PSG)$condition_imputed)

## and the areas by season for the methods part
## write.csv(table(season=sample_data(PSG)$season, area=sample_data(PSG)$area),
##           file="tables/seasonArea.csv")
## ## this is now only in-text
table(season=sample_data(PSG)$season, area=sample_data(PSG)$area)

## and save PSG for further use 
saveRDS(PSG, file="intermediate_data/PhyloSeqGenus.Rds")
