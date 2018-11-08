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

doMultiAmp <- FALSE

doTax <- FALSE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Primer test fox#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Metabarcoding/Primer_Test_Fox" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("\\d+-(Chip\\d-\\w\\d)_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[2]])

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Victors_playground/Metabarcoding/filtered/"

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
                  truncLen=c(200,200), minLen=c(200,200), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file 

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/primer_file_foxes.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primers <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
    MA <- MultiAmplicon(primers, files)
    filedir <- "/SAN/Victors_playground/Metabarcoding/stratified_files"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, n=1e+05, filedir=filedir)

    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 12)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 12)

    MA <- derepMulti(MA, mc.cores=12) 
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=12)
    MA <- mergeMulti(MA, mc.cores=12) 

    propMerged <- MultiAmplicon::calcPropMerged(MA)
        
    MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=12) 

    MA <- makeSequenceTableMulti(MA, mc.cores=12) ## FIXME in package!!!

    MA <- removeChimeraMulti(MA, mc.cores=12)

    saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/MA_final.RDS")
} else{
    MA <- readRDS("/SAN/Victors_playground/Metabarcoding/MA_final.RDS")
}

tracking <- getPipelineSummary(MA) 
## doesn't work for now

## plotPipelineSummary(tracking) 
## plotPipelineSummary(tracking) + scale_y_log10()

###Extract sequences to do taxonomic assignment 

STNC <- getSequenceTableNoChime(MA)

sequences <- unlist(lapply(STNC, colnames))
names(sequences) <- paste0("asv_", 1:length(sequences))

if(doTax){
    library(taxonomizr)
    library(taxize)

    Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
                                 "/SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta")

    ## clusters <- plotAmpliconNumbers(MA) 

###BLAST
## blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -out /SAN/Victors_playground/Metabarcoding/asv_vs_nt_fox.asn

## blast_formatter -archive /SAN/Victors_playground/Metabarcoding/asv_vs_nt_fox.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > /SAN/Victors_playground/Metabarcoding/asv_nt_fox.blttax

###Read blast result 
## we read that ouput into R blast <-
    blast <- read.csv("/SAN/Victors_playground/Metabarcoding/asv_nt_fox.blttax", header=FALSE)

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

    blt <- blt[,.(bitsum=sum(bitscore),
                  superkingdom, phylum, class, order, family, genus, species),
               by=c("query", "subject")]

    blt <- unique(blt)
    
    blt <- blt[,.(bitdiff= bitsum - max(bitsum),
                  superkingdom, phylum, class, order, family, genus, species),
               by=c("query")]

    get.unique.or.na <- function (x){
        ## unique taxa at that level excluding potential NA's 
        ux <- unique(as.character(x[!is.na(x)]))
        ## but return NA if they are not unique
        if(length(ux)==1){return(ux)} else {as.character(NA)}
    }

    genus <- blt[bitdiff>-10, .(genus=get.unique.or.na(genus)),
                 by=query]

    family <- blt[bitdiff>-20, .(family=get.unique.or.na(family)),
                  by=query]

    order <- blt[bitdiff>-30, .(order=get.unique.or.na(order)),
                 by=query]

    class <- blt[bitdiff>-40, .(class=get.unique.or.na(class)),
                 by=query]

    phylum <- blt[bitdiff>-50, .(phylum=get.unique.or.na(phylum)),
                  by=query]

    superkingdom <- blt[bitdiff>-100, .(superkingdom=get.unique.or.na(superkingdom)),
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
    seqnametab[sequences%in%dupseq,]

    seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]

    annot.list <- lapply(STNC, function (x) {
        setkey(seqnametab, sequences)
        seqnametab[colnames(x),
                   c("superkingdom", "phylum", "class", "order", "family", "genus")]
    })

    saveRDS(annot.list, file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
} else{
    annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
}


## ## Not needed anymore
## keep <- unlist(lapply(annot.list, nrow))>0
## annot.list <- annot.list[keep]
## STNC <- STNC[keep]


## name the annotation lists to have the names of the taxa 
annot.list <- lapply(seq_along(annot.list), function (i){
    an <- as.matrix(annot.list[[i]])
    rownames(an) <- colnames(STNC[[i]])
    an
})

names(annot.list) <- names(STNC)

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
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Nematoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Apicomplexa"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))


### all.annot <- Reduce(rbind, annot.list)

library(phyloseq)

## now we can add the sample information
sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/Metabarcoding_Chip1_Chip2_CS_20180824.csv")

## including the urbanization indices
urban <- read.table("/SAN/Victors_playground/Metabarcoding/Fox_urbanization_20181029.csv", sep=";", header=TRUE)
urban$X <- NULL

## But what is going on with the urbanization table? Are all these samples raccoons
## But we didn't have 47 of them, right?
table(sample.data$Sample.ID%in%urban$IZWID)

## When merging I see that all fox labels with letters can't be merged. Did you mess thes up ?
## I'll remove this letters now for merging. Please confirm that this is ther right approach!
sample.data$IZWID <- sub("[a-z]$", "", sample.data$Sample.ID)
sample.data$IZWID[sample.data$species%in%"raccoon"] <- NA

table(sample.data$IZWID%in%urban$IZWID)

## there are still three red foxes withoug IZWID!!! 
sample.data[!sample.data$IZWID%in%urban$IZWID,]
## 609, 597 and 99 missing!

## Please check! We'll move on for now!
sample.data <- merge(sample.data, urban, by="IZWID", all.x=TRUE)


## now we can match with the sequencing sample names
rownames(sample.data) <- paste0("Chip", sample.data$chip, "-", sample.data$position.chip)

## check 
table(rownames(sample.data)%in%samples)
## good! All our sequencing data can be associated with an ID in the
## sample.data table (DNA measurments etc).


## We throw out the empty amplicons only here
keep <- unlist(lapply(annot.list, nrow))>0

PS.l <- lapply(seq_along(STNC)[keep], function(i){
  phyloseq(otu_table(STNC[[i]], taxa_are_rows=FALSE),
           sample_data(sample.data[rownames(STNC[[i]]), ]),
           tax_table(annot.list[[i]]))
})


sumSeqByTax <- function (Phy, tax) {
  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
  counts$asvCount <- as.numeric(as.character(counts$asvCount))
  tapply(counts$asvCount, counts[, tax], sum)
}

readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
names(readNumByPhylum) <- names(STNC)[keep]


readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByGenus) <- names(STNC)[keep]

## ## Don't try to write out lists of unequal component length into
## ## csvs... can't work. Just continue to work happiely in R
## ## write.csv(as.data.table(names(readNumByGenus)), "/SAN/Victors_playground/Metabarcoding/genera_list_fox.csv")


################# ## HOW TO GO ON FROM HERE (Madeleine) ## ######################
## PS.l contains (will contain) all your data. It's a list of Phyloseq
## objects you can select a Phyloseq object for a single primer pair
## from this using PS.l[[n]] (where n is an integer). You can use
## lapply to apply functions over the all objects in the list.

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


## We can also merge into a single Phyloseq object forgetting the
## information which primer a taxon / ASV is from

fill <- fillSampleTables(MA)
MA@sequenceTableFilled <- fill@sequenceTableFilled

## Analyse all at once for now
ALL <- Reduce(cbind, fill@sequenceTableFilled[keep])

## Problem: over all amplicons some ASVs are identical...
table(duplicated(colnames(ALL)))

## sum up same reads over amplicons
ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## same for tax
all.tax <- Reduce(rbind, annot.list[rownames(MA)[keep]])
all.tax <- all.tax[rownames(ALL.u), ]


PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
               sample_data(sample.data[rownames(ALL), ]),
               tax_table(all.tax))

################# ## HOW TO GO ON FROM HERE (Madeleine) ## ######################
#### PS is now a single Phyloseq object over all amplicons. 

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/PhyloSeqList.Rds")
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/PhyloSeqCombi.Rds")
