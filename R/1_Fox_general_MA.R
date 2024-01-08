## Please un comment the first time you run this and re-install packages

## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)

library(MultiAmplicon)
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
library(GGally)

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

    ## and very very harsh chimera removal
    isCruelBimera <- dada2::isBimeraDenovoTable(STFU, multithread=TRUE,
                                                minSampleFraction=0.5,
                                                allowOneOff=TRUE, maxShift = 32,
                                                ignoreNNegatives=4)

    ## the same pooled over all samples
    isPooledBimera <- dada2::isBimeraDenovo(STFU, multithread=TRUE, 
                                            allowOneOff=TRUE, maxShift = 32)

    saveRDS(isCruelBimera, "/SAN/Metabarcoding/AA_Fox/cruelBimera.Rds")
    saveRDS(isPooledBimera, "/SAN/Metabarcoding/AA_Fox/pooledBimera.Rds")

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

sumPheatmap <- plotAmpliconNumbers(MAF) ### 

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

plotAmpliconNumbers(MAFinal)

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

## stargazer(modQuantExtr, modQuantSEQ, modQuantANA, type="html",
## out="tables/suppl/quant.html")
## stargazer(modQual1Extr, modQual1SEQ, modQual1ANA, type="html",
## out="tables/suppl/qual1.html")
## stargazer(modQual2Extr, modQual2SEQ, modQual2ANA, type="html",
## out="tables/suppl/qual2.html")


PS@sam_data$IDb <- rownames(PS@sam_data)
PS@sam_data$ID <- gsub("b", "", rownames(PS@sam_data))

Sdat <- as.data.frame(as(sample_data(PS), "matrix"))

DNAData <- subset(DNA, DNA$analysed)[, c("ID", "DNAng.ul", "DNA260.230", "DNA260.280")]

## overall sequencing depth
nSeq <- data.frame(nSeq=rowSums(otu_table(PS)))

Sdat <- merge(Sdat, nSeq, by.x="IDb", by.y=0)

Sdat <- merge(Sdat, DNAData, by = "ID", all=TRUE)
rownames(Sdat) <- Sdat$IDb

## they are still alinged
all(rownames(Sdat) == rownames(PS@sam_data))
PS@sam_data <- sample_data(Sdat)

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

## writing a csv export for further use in Dep. 6 and in
## collaborations. I'd always recommend to use the "phyloseq" object
## intermediate_data/PhyloSeqGenus.Rds for more control of the
## sequencing raw data though.

MoltenPSG <- psmelt(PSG)

MoltenPSG$OTU <- NULL ## removed to save space (otherwise file too large for git)
MoltenPSG$ASVn <- NULL ## same as OTU
MoltenPSG$readsF  <-  NULL
MoltenPSG$readsR  <-  NULL

write.csv(MoltenPSG, "intermediate_data/FoxPhyloSeqGenus.csv", row.names=FALSE)

## WE EXCLUDE the JUVENILES
## WE EXCLUDE everything with NAs
## remove the juveniles
PSG <- subset_samples(PSG, age %in% "adult" &
                           !is.na(sample_data(PSG)[, "tree_cover_1000m"]))

## the categories of gut content
table(tax_table(PSG)[, "category"])
### -> When we'll pick this up for future work!

## extract the Helminths
PSGHelm <- phyloseq::subset_taxa(PSG, category%in%"Helminth")

## the categories of gut content
table(tax_table(PSGHelm)[, "category"])

## now we merge Eucoleus and Capillaria
PSGHelm <- merge_taxa(PSGHelm,
                      eqtaxa = rownames(tax_table(PSGHelm))[
                          tax_table(PSGHelm)[, "genus"] %in%
                          c("Eucoleus", "Capillaria")],
                      archetype =
                          rownames(tax_table(PSGHelm))[
                              tax_table(PSGHelm)[, "genus"] %in%
                              "Eucoleus"])


tax_table(PSGHelm)[is.na(tax_table(PSGHelm))[, "genus"] ,"genus"] <- "Eucoleus"


## and save PSG for further use 
saveRDS(PSGHelm, file="intermediate_data/PhyloSeqGenus.Rds")

### Some summary data
## this regular seasons would mean only 5 summer samples... 
table(ifelse(month(sample_data(PSGHelm)$date)>11, "winter",
      ifelse(month(sample_data(PSGHelm)$date)>8, "autumn",
      ifelse(month(sample_data(PSGHelm)$date)>5, "summer",
      ifelse(month(sample_data(PSGHelm)$date)>2, "spring", "winter")))))


## and the areas by season for the methods part
table(season=sample_data(PSGHelm)$season, area=sample_data(PSGHelm)$area)

table(sample_data(PSGHelm)$sex)



### now we need to know how the environmental data for the foxes is
### correlated

sample_data(PSGHelm) %>% unclass() %>% as.data.frame() %>%
    dplyr::select(weight_kg, sex, season, condition,
                  area, tree_cover_1000m, imperv_1000m, human_fpi_1000m,
                   nSeq, DNAng.ul, DNA260.230, DNA260.280) %>%
    mutate_at(c("weight_kg",
                "tree_cover_1000m", "imperv_1000m", "human_fpi_1000m",
                "nSeq", "DNAng.ul", "DNA260.230", "DNA260.280"),
              as.numeric) %>%
    mutate_at(c("sex", "season", "condition", "area"),
              as.factor)  -> D

D %>% dplyr::select(tree_cover_1000m, imperv_1000m, human_fpi_1000m, area) %>%
    ggpairs() -> Dcor
ggsave("figures/suppl/CorrelatPedictors.png", Dcor,
        width = 12.5, height = 10.5, units = "in")

## The sampling depth == sequencing thoughput is independent of any
## technical or biological variables.

D %>%
    dplyr::select(where(is.numeric)) %>%
    gather(-nSeq, key = "var", value = "value") %>%
    ggplot(aes(x = value, y = nSeq)) +
    facet_wrap(~ var, scales = "free", ncol=3) +
    geom_point() +
    scale_y_continuous("number of sequencing reads for sample") +
    stat_smooth() -> seqPlotNum


D %>%
    dplyr::select(where(is.factor), nSeq) %>%
    gather(-nSeq, key = "var", value = "value") %>%
    ggplot(aes(x = value, y = nSeq)) +
    facet_wrap(~ var, scales = "free", ncol=2) +
    geom_boxplot(outlier.size=0) + 
    geom_jitter(width=0.3) +
    scale_y_continuous("") -> seqPlotFac


NumberSeqVarPlot <- cowplot::plot_grid(seqPlotNum, seqPlotFac, nrow = 1, ncol = 2, 
                           rel_width = c(0.8, 0.2))


ggsave("figures/suppl/NumberSeqVar.png", NumberSeqVarPlot,
        width = 25, height = 10, units = "in")


D %>%
    ggplot(aes(x = season)) +
    geom_bar() + 
    facet_wrap(~ condition) -> conditionSeason


D %>%
    ggplot(aes(x = area)) +
    geom_bar() + 
    facet_wrap(~ condition) -> conditionArea


D %>%
    ggplot(aes(x = condition , y = weight_kg)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter()  -> conditionWeight

conditionPlotFac <- cowplot::plot_grid(conditionSeason, conditionArea,
                                       conditionWeight,
                                       nrow=3, rel_heights=c(0.6, 0.6, 1))

ggsave("figures/suppl/conditionVar.png", conditionPlotFac,
       width = 8, height = 16, units = "in")


### Plot the mapping area (Figure 1)
source("R/0_Extract_Einvir_Covariates.R")
## we need fox_evncov and imperv from this first script, won't dabble
## with the details now
fox_envcov$FoxID <- paste("Fox", fox_envcov$IZW_ID)
fox_envcov <- fox_envcov[fox_envcov$FoxID%in%sample_data(PSGHelm)$ID,]


sample_data(PSGHelm)$IZW_ID[!(sample_data(PSGHelm)$IZW_ID%in%fox_envcov$IZW_ID)]

nrow(fox_envcov) ## == 141

## double check values make sense in a map
## make env cov spatial
fox_envcov_sf <- st_as_sf(fox_envcov, coords = c("coords.x1", "coords.x2"), crs = 3035)

## static map with ggplot2 + sf 
## External circle represents the values at the 1000m buffer
b <- as(extent(4400000, 4700000, 3100000, 3400000), 'SpatialPolygons')
crs(b) <- crs(imperv)
human_fpi_crop <-  st_as_stars(crop(human_fpi, b))
human_fpi_agg_1000m <- st_as_stars(terra::aggregate(crop(human_fpi, b),
                                                    fact = 10, fun = "mean"))

## shape of federal states
boundaries <- 
  st_read("input_data/VG250_Bundeslaender_esri.geojson") %>% 
  st_transform(crs = st_crs(fox_envcov_sf)) %>% 
  filter(GEN %in% c("Berlin", "Brandenburg"))


## map study area
map_study_base <- 
  ggplot(fox_envcov_sf) +
  geom_stars(data = human_fpi_crop) +
  ## state boundaries
  geom_sf(data = boundaries, fill = NA, color = "black") +
  ## 1000m buffer
  geom_sf(size = 3, shape = 21, stroke = 1.2, fill = "white", color = "white") +
  geom_sf(size = 3, shape = 21, stroke = 0, fill = "white") +
  geom_sf(aes(color = human_fpi_1000m), shape = 16, size = 3, alpha = .7) +
  ## ## middle, collection point
  geom_sf(size = 0.6, shape = 21, stroke = .8, fill = "white", color = "black") +
  geom_sf(size = 0.6, shape = 21, stroke = 0, aes(fill = human_fpi_1000m)) +
  coord_sf(xlim = c(4410000, 4650000), ylim = c(3150000, 3387000)) +
  scale_fill_gradient(low = "grey30", high = "grey96", guide = "none") +
  scale_color_scico(
    palette = "batlow", begin = .1,
    name = "Human Footprint Index (2009)", limits = c(0, 50), breaks = 1:9*5, 
    guide = guide_colorsteps(barwidth = unit(18, "lines"), barheight = unit(.6, "lines"),
                             title.position = "top", show.limits = TRUE,
                             title.hjust = 0,
                             label = FALSE)) +
    labs(x = NULL, y = NULL) +
    theme_map()

col_legend <- cowplot::get_legend(map_study_base)

fill_plot_tmp <- map_study_base +
    scale_color_scico(guide="none") + 
    scale_fill_gradient(low = "grey30", high = "grey96",
                        limits = c(0, 50), breaks = 1:9*5, 
                        guide = guide_colorsteps(barwidth = unit(18, "lines"),
                                                 barheight = unit(.6, "lines"),
                                                 title = NULL,
                                                 title.hjust = 0, show.limits = TRUE))


fill_legend <- cowplot::get_legend(fill_plot_tmp)


map_study <- map_study_base + theme(legend.position = "none") +
  ggspatial::annotation_scale(
    location = "bl", text_family = "Open Sans", text_cex = 1.2
  ) +
  ggspatial::annotation_north_arrow(location = "tr")

## map Berlin
map_berlin <- map_study_base +
  coord_sf(xlim = c(4531042, 4576603), ylim = c(3253866, 3290780)) +
  theme_void() + 
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = .8))

## overview map
sf_world <- 
  st_as_sf(rworldmap::getMap(resolution = "low")) %>% 
  st_transform(crs = st_crs(fox_envcov_sf)) %>% 
  st_buffer(dist = 0) %>% 
  dplyr::select(ISO_A2, SOVEREIGNT, LON, continent) %>% 
  mutate(area = st_area(.))

map_europe <- 
  ggplot(sf_world) +
  geom_sf(fill = "grey80", color = "grey96", lwd = .1) +
  geom_rect(
    xmin = 4430000, xmax = 4640000, ymin = 3160000, ymax = 3385000,
    color = "#212121", fill = "#a4cbb6", size = .7
  ) +
  geom_sf_text(
    data = filter(sf_world, ISO_A2 %in% c(
      "DE", "SE", "FR", "PL", "CZ", "IT", "ES", "AT", "CH", "GB", "PT", "NL", "BE", "IR", "IS"
    )),
    aes(label = ISO_A2),
    family = "Open Sans", color = "grey40", fontface = "bold", size = 4.5,
    nudge_x = 20000, nudge_y = -10000
  ) +
  ggspatial::annotation_scale(
    location = 'tr', text_family = "Open Sans", text_cex = 1.2
  ) +
  coord_sf(xlim = c(2650000, 5150000), ylim = c(1650000, 5100000)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-10, 30, by = 10)) +
  labs(x = NULL, y = NULL) +
  theme_map() +
  theme(panel.ontop = FALSE,
        panel.grid.major = element_line(color = "grey75", linetype = "15", linewidth = .3))

map_globe <- d6berlin::globe(col_earth = "grey80", col_water = "grey96", bg = TRUE)

### combined map
map_overview <- map_europe +
    labs(tag = "a") +
    inset_element(map_globe, .02, .75, .59, 1, align_to = "plot")

map_foo <- map_study +
    inset_element(map_berlin + ggtitle("Berlin"), .14, .1, .5, 0.48,
                  align_to = "plot")

final_legend <- cowplot::plot_grid(NULL, col_legend, NULL, 
                                   NULL, fill_legend, NULL,
                                   rel_widths = c(0.13, 0.87, 0.13),
                                   nrow=2)

map_bar <- cowplot::plot_grid(final_legend,  map_foo, 
                              ncol = 1, rel_heights = c(0.15, 1.2)) +
    labs(tags = "b")

m <- cowplot::plot_grid(map_overview, map_bar,
                        ncol = 2, rel_widths = c(0.78, 1),
                        align = "h", axis = "t")

ggsave("figures/map_study_overview_multi.png", m,
       width = 11.5, height = 7, bg = "white", dpi = 600)
