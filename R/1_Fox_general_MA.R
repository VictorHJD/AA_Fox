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

sample.data$IZW_ID <- as.vector(sample.data$IZW_ID)

### We have to fix the IZW sample names!
setdiff(sample_names(PS), sample.data$IZW_ID)
## they have these additional "b"s

## just add "b"s to all those
sample.data$IZW_ID[!sample.data$IZW_ID%in%sample_names(PS)] <-
    paste0(sample.data$IZW_ID[!sample.data$IZW_ID%in%sample_names(PS)], "b")

rownames(sample.data) <- sample.data$IZW_ID

## align and cbind to get the combinded sample data
PS@sam_data <- sample_data(cbind(PS@sam_data, sample.data[rownames(sample_data(PS)), ]))

## make weight numeric
PS@sam_data[, "weight_kg"] <- as.numeric(unlist(PS@sam_data[, "weight_kg"] ))

## now directly in the repository (for reproducibilty)
saveRDS(PS, file="intermediate_data/PhyloSeqCombi.Rds")

###For primer analysis (Victor), still stored on our server 
## saveRDS(PS.l, file="/SAN/Metabarcoding/AA_Fox/PhyloSeqList.Rds") 

## For Caro and the paper, previously on the server, 
## saveRDS(PS, file="/SAN/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds")


### Adding Categories to the taxonomy information!
### We have to addd this to genus agglommerted data!
## collapse to genus level
PSG <- phyloseq::tax_glom(PS, "genus")

## We take the helminth traits information from manually curated our
## input data
traits <- read.csv("input_data/helminth_traits.csv")

traits %>%
    column_to_rownames("t.genus")  -> traits 

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
                              "Mollusca")&
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

table(tax_table(PSG)[, "category"])
## and save PSG for further use 
saveRDS(PSG, file="intermediate_data/PhyloSeqGenus.Rds")
