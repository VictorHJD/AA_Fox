## Please uncomment the first time you run this and re-install packages

## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(MultiAmplicon)
library(ggplot2)
library(data.table)

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

##Multi amplicon pipeline
if(doMultiAmp){
  MAF <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Fox/Stratified_files_complete"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MAF <- sortAmplicons(MAF, n=1e+05, filedir=filedir)
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MAF)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MAF)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  MAF <- derepMulti(MAF, mc.cores=12) 
  MAF <- dadaMulti(MAF, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  MAF <- mergeMulti(MAF, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MAF)
  
  MAF <- mergeMulti(MAF, justConcatenate=propMerged<0.8, mc.cores=12) 
  
  MAF <- makeSequenceTableMulti(MAF, mc.cores=12) ## FIXME in package!!!
  
  MAF <- removeChimeraMulti(MAF, mc.cores=12)
  
  saveRDS(MAF, "/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS")
} else{
  MAF <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS") ###START from here now! 
}

trackingF <- getPipelineSummary(MAF) 
## doesn't work for now

plotAmpliconNumbers(MAF) ### Second batch of sequences did not sum any useful data :(


## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

###Extract sequences to do taxonomic assignment 

STNCF <- getSequenceTableNoChime(MAF)

sequences <- unlist(lapply(STNCF, colnames))
names(sequences) <- paste0("asv_", 1:length(sequences))

if(doTax){
  library(taxonomizr)
  library(taxize)
  
  Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
                              "/SAN/Victors_playground/Metabarcoding/AA_Fox/FoxRun_seq_complete.fasta")
  
  clusters <- plotAmpliconNumbers(MAF) 
  
  ###BLAST
  ## blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Metabarcoding/AA_Fox/FoxRun_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -out /SAN/Victors_playground/Metabarcoding/AA_Fox/asv_vs_nt_foxfinal.asn
  
  ## blast_formatter -archive /SAN/Victors_playground/Metabarcoding/AA_Fox/asv_vs_nt_foxfinal.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > /SAN/Victors_playground/Metabarcoding/AA_Fox/asv_nt_foxfinal.blttax
  
  ###Read blast result 
  ## we read that ouput into R blast <-
  blast <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Fox/asv_vs_nt_FoxComplete.blttax", header=FALSE)
  
  names(blast) <- c("query", "subject", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                    "bitscore", "staxid")
  blast <- as.data.table(blast)
  blast$staxid <- as.character(blast$staxid)
  
  read.nodes.sql("/SAN/db/taxonomy/nodes.dmp",
                 "/SAN/db/taxonomy/taxonomizr.sql")
  read.names.sql("/SAN/db/taxonomy/names.dmp",
                 "/SAN/db/taxonomy/taxonomizr.sql")
  
  blast.tax <- getTaxonomy(unique(blast$staxid),
                           "/SAN/db/taxonomy/taxonomizr.sql")
  
  blast.tax <- as.data.table(blast.tax, keep.rownames="staxid")
  blast.tax$staxid <- gsub("\\s*", "", blast.tax$staxid)
  
  blt <- merge(blast, blast.tax, by="staxid", all=TRUE)
  
  ## ## ## We need to be more clever if we want to use multiple
  ## ## ## hsps, this does not work for whole genome subjects eg.
  ### blt <- blt[,.(bitsum=sum(bitscore),
  ###              superkingdom, phylum, class, order, family, genus, species),
  ###           by=c("query", "subject")]
  
  ###    blt <- unique(blt)
  
  blt <- blt[,.(bitdiff= bitscore - max(bitscore),
                superkingdom, phylum, class, order, family, genus, species),
             by=c("query")]
  
  get.unique.or.na <- function (x){
    ## unique taxa at that level excluding potential NA's 
    ux <- unique(as.character(x[!is.na(x)]))
    ## but return NA if they are not unique
    if(length(ux)==1){return(ux)} else {as.character(NA)}
  }
  
  genus <- blt[bitdiff>-2, .(genus=get.unique.or.na(genus)),
               by=query]
  
  family <- blt[bitdiff>-7, .(family=get.unique.or.na(family)),
                by=query]
  
  order <- blt[bitdiff>-12, .(order=get.unique.or.na(order)),
               by=query]
  
  class <- blt[bitdiff>-20, .(class=get.unique.or.na(class)),
               by=query]
  
  phylum <- blt[bitdiff>-30, .(phylum=get.unique.or.na(phylum)),
                by=query]
  
  superkingdom <- blt[bitdiff>-50, .(superkingdom=get.unique.or.na(superkingdom)),
                      by=query]
  
  annot <- cbind(superkingdom[,c("query", "superkingdom")],
                 phylum[,"phylum"],
                 class[,"class"],
                 order[,"order"],
                 family[,"family"],
                 genus[,"genus"])
  
  seqnametab <- as.data.table(cbind(query=names(sequences), sequences))
  seqnametab <- merge(seqnametab, annot)
  
  dupseq <- seqnametab$sequences[duplicated(seqnametab$sequences)]
  
  seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]
  
  annot.list <- lapply(STNCF, function (x) {
    setkey(seqnametab, sequences)
    seqnametab[colnames(x),
               c("superkingdom", "phylum", "class", "order", "family", "genus")]
  })
  
  saveRDS(annot.list, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_blast_tax_complete.Rds")
} else{
  annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_blast_tax_complete.Rds")
}


## ## Not needed anymore
## keep <- unlist(lapply(annot.list, nrow))>0
## annot.list <- annot.list[keep]
## STNC <- STNC[keep]


## name the annotation lists to have the names of the taxa 
annot.list <- lapply(seq_along(annot.list), function (i){
  an <- as.matrix(annot.list[[i]])
  rownames(an) <- colnames(STNCF[[i]])
  an
})

names(STNCF)<-gsub(pattern = "-", replacement = "_", x= names(STNCF))
names(STNCF)<-gsub(pattern = " ", replacement = "", x= names(STNCF))

names(annot.list) <- names(STNCF)

phylalist <- lapply(annot.list, function (x) {
  if(nrow(x)>0){
    table(x[, "phylum"])
  }
})


tabulate.taxa <- function(taxtab, taxon, phylumsubset){
  if(nrow(taxtab)>0){
    t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
    if(!is.null(ncol(t))){
      table(t[, taxon])
    } else {NULL} 
  }else {NULL} 
}


## Tabulate by specific phylum
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Cestoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Nematoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Apicomplexa"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "family", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Amphibia"))

### all.annot <- Reduce(rbind, annot.list)

library(phyloseq)

## now we can add the sample information
sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_data.csv")


## including the urbanization indices
#urban <- read.table("/SAN/Victors_playground/Metabarcoding/Fox_urbanization_20181029.csv", sep=";", header=TRUE)
#urban$X <- NULL

## But what is going on with the urbanization table? Are all these samples raccoons
## But we didn't have 47 of them, right?
#table(sample.data$Sample.ID%in%urban$IZWID)

## When merging I see that all fox labels with letters can't be merged. Did you mess thes up ?
## I'll remove this letters now for merging. Please confirm that this is ther right approach!
#sample.data$IZWID <- sub("[a-z]$", "", sample.data$Sample.ID)
#sample.data$IZWID[sample.data$species%in%"raccoon"] <- NA

#table(sample.data$IZWID%in%urban$IZWID)

## there are still three red foxes withoug IZWID!!! 
#sample.data[!sample.data$IZWID%in%urban$IZWID,]
## 609, 597 and 99 missing!

## Please check! We'll move on for now!
#sample.data <- merge(sample.data, urban, by="IZWID", all.x=TRUE)


## now we can match with the sequencing sample names
#rownames(sample.data) <- paste0("Chip", sample.data$chip, "-", sample.data$position.chip)

## check 
#table(rownames(sample.data)%in%samples)
## good! All our sequencing data can be associated with an ID in the
## sample.data table (DNA measurments etc).

rownames(sample.data) <- sample.data$IZW_ID ###take sample id as rowname an make the function below work ;)  

## We throw out the empty amplicons only here
keep <- unlist(lapply(annot.list, nrow))>0

rownames(STNCF)

PS.l <- lapply(seq_along(STNCF)[keep], function(i){
  phyloseq(otu_table(STNCF[[i]], taxa_are_rows=FALSE),
           sample_data(sample.data[rownames(STNCF[[i]]), ]),
           tax_table(annot.list[[i]]))
})


sumSeqByTax <- function (Phy, tax) {
  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
  counts$asvCount <- as.numeric(as.character(counts$asvCount))
  tapply(counts$asvCount, counts[, tax], sum)
}

readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
names(readNumByPhylum) <- names(STNCF)[keep]


readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByGenus) <- names(STNCF)[keep]


readNumByFamily <- lapply(PS.l, sumSeqByTax, "family") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByFamily) <- names(STNCF)[keep]

####
fill <- fillSampleTables(MAF)
MAF@sequenceTableFilled <- fill@sequenceTableFilled


## Analyse all at once for now
ALL <- Reduce(cbind, fill@sequenceTableFilled[keep])

## Problem: over all amplicons some ASVs are identical...
table(duplicated(colnames(ALL)))

## sum up same reads over amplicons
ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## same for tax
all.tax <- Reduce(rbind, annot.list[rownames(MAF)[keep]])
all.tax <- all.tax[rownames(ALL.u), ]

PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
               sample_data(sample.data[rownames(ALL), ]),
               tax_table(all.tax))

prune_both_zero <- function (ps) {
  p <- prune_samples(sample_sums(ps) > 0 , ps)
  prune_taxa(taxa_sums(p) > 0 , p)
}

PS <- prune_both_zero(PS)
PS.l <- lapply(PS.l, prune_both_zero)

################# ## HOW TO GO ON FROM HERE ## ######################
#### PS is now a single Phyloseq object over all amplicons. 

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqList.Rds") ###For primer analysis (Victor)
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqCombi.Rds") ###For Fox analysis (Caro and Sophia)
