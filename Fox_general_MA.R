## Please uncomment the first time you run this and re-install packages

## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
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

doMultiAmp <- TRUE

doTax <- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run foxes#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/2018_22_fox_all/" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

#samples <- gsub("s\\d+_", "\\1", basename(fastqF))
samples <- gsub("_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
#samples<- gsub("-", "_", basename(samples))

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
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Fox/primer_file_foxes.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)
MAF <- MultiAmplicon(primer, files)

##Multi amplicon pipeline
if(doMultiAmp){
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/Stratified_files_new"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MAF <- sortAmplicons(MAF, n=1e+05, filedir=filedir)
  ##Old pipeline
  #filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/Stratified_files_complete"
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MAF)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MAF)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  #MAF <- derepMulti(MAF, mc.cores=12) 
  MAF <- dadaMulti(MAF, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  MAF <- mergeMulti(MAF, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MAF)
  
  summary(propMerged)
  table(propMerged<0.8)
  
  MAF <- mergeMulti(MAF, justConcatenate=propMerged<0.8, mc.cores=12) 
  
  MAF <- makeSequenceTableMulti(MAF, mc.cores=12)
  
  MAF <- removeChimeraMulti(MAF, mc.cores=12)
  
  saveRDS(MAF, "/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS")
} else{
  MAF <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS") ###START from here now! 
}

## When loading an old MA object that lacks sample data, simply:
MAF <- addSampleData(MAF)

trackingF <- getPipelineSummary(MAF) 
## doesn't work for now

plotAmpliconNumbers(MAF) ### 


## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

###Extract sequences to do taxonomic assignment 

STNCF <- getSequenceTableNoChime(MAF)

sequences <- unlist(lapply(STNCF, colnames))
names(sequences) <- paste0("asv_", 1:length(sequences))

###New taxonomic assignment 
#MAF <- blastTaxAnnot(MAF,  dataBaseDir = Sys.getenv("BLASTDB"), negative_gilist = "/SAN/db/blastdb/uncultured.gi", num_threads = 15)

if (doTax){ ## simply save the blast files, that's even faster than
  ## setting doTax to FALSE and re-loading the object
  MAF2 <- blastTaxAnnot(MAF,  
                        negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                        db = "/SAN/db/blastdb/nt/nt",
                        infasta = "/SAN/Victors_playground/Metabarcoding/AA_Fox/in.fasta",
                        outblast = "/SAN/Victors_playground/Metabarcoding/AA_Fox/out.blt",
                        num_threads = 22)
  saveRDS(MAF2, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF2.Rds")
} else {
  MAF2 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF2.Rds")
}

##Add real sample data
sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_data.csv", dec=",", stringsAsFactors=FALSE)
sample.data$IZW_ID <- as.vector(sample.data$IZW_ID)
rownames(sample.data) <- sample.data$IZW_ID
MAF3 <- addSampleData(MAF2, sample.data)
saveRDS(MAF3, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF3.Rds")

MAF3 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF3.Rds")

###Couple of checks before phyloseq
lapply(getTaxonTable(MAF3), function (x) table(as.vector(x[, "phylum"])))
lapply(getTaxonTable(MAF3), function (x) table(as.vector(x[, "genus"])))
lapply(getTaxonTable(MAF3), function (x) table(as.vector(x[, "species"])))

##to phyloseq

PS.l <- toPhyloseq(MAF3, samples=colnames(MAF3), multi2Single=FALSE)

PS <- toPhyloseq(MAF3, samples=colnames(MAF3), multi2Single=TRUE)

#pdf(file = "~/AA_Primer_evaluation/Figures/Fox_Rowreads.pdf", width = 10, height = 20)
plotAmpliconNumbers(MAF3, cluster_cols= T, cluster_row=F,cutree_cols= 2)
#dev.off()

saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqList.Rds") ###For primer analysis (Victor)
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds") ###For Fox analysis (Caro and Sophia)
