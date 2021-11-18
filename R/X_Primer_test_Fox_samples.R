## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

## library(MultiAmplicon)

## using the development version
devtools::load_all("../MultiAmplicon")

library(ggplot2)
library(reshape)
library(taxize)
library(parallel)
library(phyloseq)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE

doMultiAmp <- TRUE

doTax <- TRUE
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
    filedir <- "/SAN/Metabarcoding/AA_Fox/stratified_files"
    if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
    MA <- sortAmplicons(MA, filedir=filedir)

    errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                         verbose=0, multithread = 48)
    errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                        verbose=0, multithread = 48)

    MA <- derepMulti(MA, mc.cores=48) 
    MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                    verbose=0, mc.cores=48)
    MA <- mergeMulti(MA, mc.cores=48) 

    propMerged <- MultiAmplicon::calcPropMerged(MA)
        
    MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=48) 

    MA <- makeSequenceTableMulti(MA, mc.cores=12) ## FIXME in package!!!

    MA <- removeChimeraMulti(MA, mc.cores=12)

    saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/MAFR_Test.RDS")
} else{
    MA <- readRDS("/SAN/Victors_playground/Metabarcoding/MAFR_Test.RDS")
}

plotAmpliconNumbers(MA, cluster_cols= T)
tracking <- getPipelineSummary(MA) 
## doesn't work for now

###Taxonomic annotation
MA2 <- blastTaxAnnot(MA,  
                      negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                      db = "/SAN/db/blastdb/nt/nt",
                      infasta = "/SAN/Victors_playground/Metabarcoding/in.fasta",
                      outblast = "/SAN/Victors_playground/Metabarcoding/out.blt",
                      num_threads = 22)
saveRDS(MA2, file="/SAN/Victors_playground/Metabarcoding/MAR2.Rds")

MA2 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/MAR2.Rds")

lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "phylum"])))
lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "genus"])))
lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "species"])))

PS.l <- toPhyloseq(MA2, samples=colnames(MA2), multi2Single=FALSE)

PS <- toPhyloseq(MA2, samples=colnames(MA2), multi2Single=TRUE)

taxa_sums(MA2) 


## plotPipelineSummary(tracking) 
## plotPipelineSummary(tracking) + scale_y_log10()

###Extract sequences to do taxonomic assignment 

#STNC <- getSequenceTableNoChime(MA)

#sequences <- unlist(lapply(STNC, colnames))
#names(sequences) <- paste0("asv_", 1:length(sequences))

#if(doTax){
#    library(taxonomizr)
#    library(taxize)

    ## Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
    ##                             "/SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta")

    ## clusters <- plotAmpliconNumbers(MA) 

###BLAST
## blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -out /SAN/Victors_playground/Metabarcoding/asv_vs_nt_fox.asn

## blast_formatter -archive /SAN/Victors_playground/Metabarcoding/asv_vs_nt_fox.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > /SAN/Victors_playground/Metabarcoding/asv_nt_fox.blttax

###Read blast result 
## we read that ouput into R blast <-
#    blast <- read.csv("/SAN/Victors_playground/Metabarcoding/asv_nt_fox.blttax", header=FALSE)

#    names(blast) <- c("query", "subject", "pident", "length", "mismatch",
  #                    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
  #                    "bitscore", "staxid")
#    blast <- as.data.table(blast)
 #   blast$staxid <- as.character(blast$staxid)
    
 #   read.nodes.sql("/SAN/db/taxonomy/nodes.dmp",
 #                  "/SAN/db/taxonomy/taxonomizr.sql")
#    read.names.sql("/SAN/db/taxonomy/names.dmp",
#                   "/SAN/db/taxonomy/taxonomizr.sql")

#    blast.tax <- getTaxonomy(unique(blast$staxid),
#                             "/SAN/db/taxonomy/taxonomizr.sql")

 #   blast.tax <- as.data.table(blast.tax, keep.rownames="staxid")
 #   blast.tax$staxid <- gsub("\\s*", "", blast.tax$staxid)
    
 #   blt <- merge(blast, blast.tax, by="staxid", all=TRUE)

    ## ## ## We need to be more clever if we want to use multiple
    ## ## ## hsps, this does not work for whole genome subjects eg.
    ### blt <- blt[,.(bitsum=sum(bitscore),
    ###              superkingdom, phylum, class, order, family, genus, species),
    ###           by=c("query", "subject")]

    ###    blt <- unique(blt)
    
#    blt <- blt[,.(bitdiff= bitscore - max(bitscore),
#                  superkingdom, phylum, class, order, family, genus, species),
 #              by=c("query")]

 #   get.unique.or.na <- function (x){
        ## unique taxa at that level excluding potential NA's 
 #       ux <- unique(as.character(x[!is.na(x)]))
        ## but return NA if they are not unique
#        if(length(ux)==1){return(ux)} else {as.character(NA)}
#    }

#    genus <- blt[bitdiff>-2, .(genus=get.unique.or.na(genus)),
 #                by=query]

#    family <- blt[bitdiff>-7, .(family=get.unique.or.na(family)),
#                  by=query]

#    order <- blt[bitdiff>-12, .(order=get.unique.or.na(order)),
#                 by=query]

#    class <- blt[bitdiff>-20, .(class=get.unique.or.na(class)),
#                 by=query]

#    phylum <- blt[bitdiff>-30, .(phylum=get.unique.or.na(phylum)),
#                  by=query]

#    superkingdom <- blt[bitdiff>-50, .(superkingdom=get.unique.or.na(superkingdom)),
#                        by=query]

#    annot <- cbind(superkingdom[,c("query", "superkingdom")],
#                   phylum[,"phylum"],
#                   class[,"class"],
#                   order[,"order"],
#                   family[,"family"],
#                   genus[,"genus"])

#    seqnametab <- as.data.table(cbind(query=names(sequences), sequences))
#    seqnametab <- merge(seqnametab, annot)

#    dupseq <- seqnametab$sequences[duplicated(seqnametab$sequences)]
    
#    seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]

#    annot.list <- lapply(STNC, function (x) {
#        setkey(seqnametab, sequences)
#        seqnametab[colnames(x),
#                   c("superkingdom", "phylum", "class", "order", "family", "genus")]
#    })

#    saveRDS(annot.list, file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
#} else{
#    annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
#}


## ## Not needed anymore
## keep <- unlist(lapply(annot.list, nrow))>0
## annot.list <- annot.list[keep]
## STNC <- STNC[keep]


## name the annotation lists to have the names of the taxa 
#annot.list <- lapply(seq_along(annot.list), function (i){
#    an <- as.matrix(annot.list[[i]])
#    rownames(an) <- colnames(STNC[[i]])
#    an
#})

#names(annot.list) <- names(STNC)

#phylalist <- lapply(annot.list, function (x) {
#     if(nrow(x)>0){
#         table(x[, "phylum"])
#     }
#})


#tabulate.taxa <- function(taxtab, taxon, phylumsubset){
#    if(nrow(taxtab)>0){
#        t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
#        if(!is.null(ncol(t))){
#            table(t[, taxon])
#        } else {NULL} 
#    }else {NULL} 
#}


## Tabulate by specific phylum
#lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Chordata"))
#lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Nematoda"))
#lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Apicomplexa"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))

#lapply(annot.list, function (x) tabulate.taxa(x, "family", "Chordata"))
#lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))

### all.annot <- Reduce(rbind, annot.list)



## now we can add the sample information
sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/Metabarcoding_Chip1_Chip2_CS_20180824.csv")

sample.dataFox <- read.csv("/SAN/Victors_playground/Metabarcoding/Fox_Info.csv")

sample.data <- merge(sample.data, sample.dataFox, all.x=TRUE)

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
#keep <- unlist(lapply(annot.list, nrow))>0

#PS.l <- lapply(seq_along(STNC)[keep], function(i){
 # phyloseq(otu_table(STNC[[i]], taxa_are_rows=FALSE),
 #          sample_data(sample.data[rownames(STNC[[i]]), ]),
 #          tax_table(annot.list[[i]]))
#})


#sumSeqByTax <- function (Phy, tax) {
 # counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
#  counts$asvCount <- as.numeric(as.character(counts$asvCount))
#  tapply(counts$asvCount, counts[, tax], sum)
#}

#readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
#names(readNumByPhylum) <- names(STNC)[keep]


#readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
#names(readNumByGenus) <- names(STNC)[keep]

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

#fill <- fillSampleTables(MA)
#MA@sequenceTableFilled <- fill@sequenceTableFilled

## Analyse all at once for now
#ALL <- Reduce(cbind, fill@sequenceTableFilled[keep])



## Problem: over all amplicons some ASVs are identical...
#table(duplicated(colnames(ALL)))

## sum up same reads over amplicons
#ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## same for tax
#all.tax <- Reduce(rbind, annot.list[rownames(MA)[keep]])
#all.tax <- all.tax[rownames(ALL.u), ]


#PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
#               sample_data(sample.data[rownames(ALL), ]),
#               tax_table(all.tax))

#prune_both_zero <- function (ps) {
#    p <- prune_samples(sample_sums(ps) > 0 , ps)
#    prune_taxa(taxa_sums(p) > 0 , p)
#}

#PS <- prune_both_zero(PS)
#PS.l <- lapply(PS.l, prune_both_zero)

################# ## HOW TO GO ON FROM HERE (Madeleine) ## ######################
#### PS is now a single Phyloseq object over all amplicons. 

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/PhyloSeqList.Rds")
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/PhyloSeqCombi.Rds")

#######Primer evaluation Victor#####
PS.l <- readRDS("/SAN/Victors_playground/Metabarcoding/PhyloSeqList.Rds")

MA2@taxonTable

readNumByPhylum<- lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "phylum"])))

rowSums(getRawCounts(MA)) ### Total amount of read per primer pair 
rawcounts <- rowSums(getRawCounts(MA))
rawcounts[order(rawcounts)]

rawcounts <- data.frame(rawcounts)
rawcounts[,2] <- rownames(rawcounts)
colnames(rawcounts) <- c("Raw_counts", "Primer_name")
rownames(rawcounts) <- c(1:nrow(rawcounts))
rawcounts <- data.frame(Primer_name = rawcounts$Primer_name, Raw_counts = rawcounts$Raw_counts) ###change the order of the columns
rawcounts$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts$Primer_name) 
rawcounts$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts$Primer_name) 

colSums(getRawCounts(MA)) ## Total amount of reads per sample 

sum(rawcounts$Raw_counts)

readNumByPhylum <- lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "phylum"])))
readNumByGenus <- lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "genus"])))
readNumByFamily <- lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "family"])))
readNumByOrder <- lapply(getTaxonTable(MA2), function (x) table(as.vector(x[, "order"])))

#######Nice data frame Phylum##### 
library(dplyr)

AbPhy <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: fro every element in the list ...
{ 
  phyla <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    phyla[1,1] <- 0    ### Add a zero in the first column
    phyla[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    phyla <- as.data.frame((readNumByPhylum[[i]]))  ###Make a data frame with the data included in each element of the list 
    phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
  }
  
  phyla[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhy <- rbind(AbPhy, phyla) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_name = AbPhy$Primer_name, Phyla = AbPhy$Phyla, Reads = AbPhy$Reads) ###change the order of the columns

AbPhy %>%
  group_by(Primer_name) %>% 
  mutate(Total_reads = sum(Reads)) -> AbPhy ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
 
Relative_abundance = AbPhy$Reads/AbPhy$Total_reads ### create a vector with the result of the operation 

AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

#unique(AbPhy$Primer_name)

pinfo <- read.csv("/SAN/Victors_playground/Metabarcoding/Primer_information.csv") 
pinfo <- pinfo[,c("Primer_pair_name", "Target_gene", "Group", "Amplicon_expected", "Access_array_Status")]

AbPhy$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_name) ### check if the primer names have extra spaces

#setdiff(AbPhy$Primer_name, pinfo$Primer_pair_name)
#setdiff(pinfo$Primer_pair_name, AbPhy$Primer_name)

AbPhy <- merge(AbPhy, pinfo, by.x= "Primer_name", by.y = "Primer_pair_name") ###merge the selected information with the origial data frame created 

write.csv(AbPhy, file = "/SAN/Victors_playground/Metabarcoding/Phyla_list_experiment_0618_Fox.csv")

#############Nice data frame Genus##### 
library(dplyr)

AbGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Reads", "Genus", "Primer_name") ### change the names for the columns 
  AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_name = AbGen$Primer_name, Genus = AbGen$Genus, Reads = AbGen$Reads) ###change the order of the columns

AbGen %>%
  group_by(Primer_name) %>% 
  mutate(Total_reads = sum(Reads)) -> AbGen ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbGen$Reads/AbGen$Total_reads ### create a vector with the result of the operation 

AbGen[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbGen)[5] <- "Relative_abundance" ### Cahnge the name of the column 

#unique(AbPhy$Primer_name)

#pinfo <- read.csv("/SAN/Victors_playground/Metabarcoding/Primer_information.csv") 
#pinfo <- pinfo[,c("Primer_pair_name", "Target_gene", "Group", "Amplicon_expected", "Access_array_Status")]

AbGen$Primer_name <- gsub(pattern = " ", replacement = "", x = AbGen$Primer_name) ### check if the primer names have extra spaces

#setdiff(AbPhy$Primer_name, pinfo$Primer_pair_name)
#setdiff(pinfo$Primer_pair_name, AbPhy$Primer_name)

AbGen <- merge(AbGen, pinfo, by.x= "Primer_name", by.y = "Primer_pair_name") ###merge the selected information with the origial data frame created 

write.csv(AbGen, file = "/SAN/Victors_playground/Metabarcoding/Genus_list_experiment_0618_Fox.csv")

#########Plots for primers#####

library("ggplot2")

##Total reads
richness <- data.frame()

for (i in 1:length(PS.l)){
  
  a <- data.frame()
  
  b <- data.frame(estimate_richness(PS.l[[i]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
  
  a<-colMeans(b)
  
  richness<- rbind(richness, a)
}

richness[,7] <- names(STNC[keep])

colnames(richness) <- c("Observed","Chao1","se_Chao1", "Shannon", "Simpson", "InvSimpson", "Primer_name")

nsamp <- data.frame()

nsamp[1:55,1] <-names(STNC[keep])

for (i in 1:length(PS.l))
  nsamp[i,2]<- nrow(sample_data(PS.l[[i]]))

colnames(nsamp) <- c("Primer_name", "Number_samples")

richness$Primer_name <- gsub(pattern = " ", replacement = "", x = richness$Primer_name) ### check if the primer names have extra spaces
nsamp$Primer_name <- gsub(pattern = " ", replacement = "", x= nsamp$Primer_name)

#primers$Primer_name <- gsub(pattern = " ", replacement = "", x= primers$Primer_name)

#setdiff(richness$Primer_name, nsamp$Primer_name)

richness <- merge(richness, nsamp, by= "Primer_name")

richness <- richness %>% distinct(Primer_name, .keep_all = TRUE)

richness$Primer_name <- gsub(pattern = "-", replacement = "_", x = richness$Primer_name)

richness <- merge(richness, rawcounts, "Primer_name", all= T) %>% distinct(Primer_name, .keep_all = TRUE)

richness <- merge(richness, primerInput, "Primer_name")

richness <- richness[which(richness$Raw_counts>8),] 

Div <- richness$Shannon
Rich <- log(richness$Observed + 1)
evenness = Div/Rich
richness$Evenness = evenness

###Plot raw counts by primer pair 
ggplot(richness, aes(x= reorder(Primer_name, -Raw_counts), y=Raw_counts, color= Gen, fill= Gen)) + 
  geom_bar(stat = "identity")+
  scale_y_log10(name = "log10 Total sequencing reads")+
  theme_minimal() +
  labs(x= "Primer name")+
  geom_hline(yintercept=100, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))


###Amplicon size vs raw counts
ggplot(richness, aes(x = Expected, y = Raw_counts, color= Target), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Expected amplicon size") +
  #scale_y_continuous(name = "Sequencing reads") + #, limits=c(0.95, 1.60)) + ggtitle("Oocysts L/W ratio by group")+
  scale_y_log10(name = "log10 Total sequencing reads")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5, aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  geom_hline(yintercept=100, linetype="dashed", color = "black") + 
  theme(text = element_text(size=20))

adjrsq <- summary(lm(pinfo$Amplicon_expected ~ pinfo$Raw_counts))$adj.r.squared



##Phylum
ggplot(data=AbPhy, aes(x= Primer_name, y= Relative_abundance, fill= Phyla)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 


###Genus
ggplot(data=AbGen, aes(x= Primer_name, y= Relative_abundance, fill= Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="down") + guides(fill=guide_legend(nrow=5)) 


length((unique(na.omit(AbPhy$Phyla))))
length((unique(na.omit(AbGen$Genus))))

###Bacterial biome
BacP16S <- subset(AbPhy, Target_gene == "16S" & Group == "Bacteria_biome")

ggplot(data = BacP16S , aes(x = Primer_name, y = Relative_abundance, fill = Phyla)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Genl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) 

###18S
Phy18S <- subset(AbPhy, AbPhy$Target_gene == "18S")
ggplot(data = Phy18S, aes(x = Primer_name, y = Relative_abundance, fill = Phyla)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Genl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) 



###Plant

unique(GenPlant$Genus)
count(subset(AbGen, AbGen$Primer_name == "trnL(UAA)g-65_F.trnL(UAA)h-65_R"))
count(subset(AbPhy, AbPhy$Primer_name == "trnL(UAA)g-65_F.trnL(UAA)h-65_R"))

subset(AbGen, AbGen$Group == "Plant_Diet")
GenPlant <- subset(AbGen, AbGen$Group == "Plant_Diet")
ggplot(data = GenPlant, aes(x = Primer_name, y = Relative_abundance, fill = Genus)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=25)) 

###Diet 16S
subset(pinfo, pinfo$Group == "Animal_Diet")
GenDiet16 <- subset(AbGen,  AbGen$Target_gene != "CytB" & AbGen$Group == "Animal_Diet")
ggplot(data = GenDiet16, aes(x = Primer_name, y = Relative_abundance, fill = Genus)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 

GenDiet <- subset(AbGen, AbGen$Group == "Animal_Diet")
ggplot(data = GenDiet, aes(x = Primer_name, y = Relative_abundance, fill = Genus)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 

###COI
count(subset(pinfo, pinfo$Target_gene == "COI"))
GenCOI <- subset(AbGen,  AbGen$Target_gene == "COI")
ggplot(data = GenCOI, aes(x = Primer_name, y = Relative_abundance, fill = Genus)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 

pinfo.coi <- subset(pinfo, pinfo$Target_gene == "COI")
ggplot(pinfo.coi, aes(x= reorder(Primer_name, -Raw_counts), y=Raw_counts, color= Target_gene, fill= Target_gene)) + 
  geom_bar(stat = "identity")+
  #scale_y_log10(name = "log10 Total sequencing reads")+
  scale_y_continuous(name = "Total sequencing reads") +
  theme_minimal() +
  labs(x= "Primer name")+
  geom_hline(yintercept=100, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))




coord_polar() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") + ylab("")


######

library("phyloseq")
library("ggplot2")
library("vegan")
library("r")

PS.l <- readRDS(file= "/SAN/Victors_playground/Metabarcoding/PhyloSeqList.Rds")
names(PS.l) <- names(STNC)[keep] 
PS.l



prune_both_zero <- function (ps) {
  p <- prune_samples(sample_sums(ps) > 0 , ps)
  prune_taxa(taxa_sums(p) > 0 , p)
}

PS <- prune_both_zero(PS)
PS.l <- lapply(PS.l, prune_both_zero)


nameOtuByTax <- function(ps, taxon="genus"){
  otable <- otu_table(ps)
  rownames(otable) <- make.unique(as.character(tax_table(ps)[, "genus"]))
  otable
}


make.unique(as.character(tax_table(PS.l[["Ave12F-91_F.Ave12R-92_R"]])[, "genus"]))

lapply(PS.l, function(x) { 
  make.unique(as.character(tax_table(x[, "genus"])))
}
)


taxa <- lapply(PS.l, function(x, taxon="genus"){
  otable <- otu_table(x)
  rownames(otable) <- make.unique(as.character(tax_table(x)[, "genus"]))
  otable
}
)

##Rarefaction curves
##made individual curves
lapply(PS.l, function(x){
  rarecurve(otu_table(x), cex= 0.5)
}
  )                     

min_samp <- lapply(PS.l, function(x){
  min(sample_sums(x))
}
) 

set.seed(123)


PS.pri_rare <- lapply(PS.l, function(x){
  rarefy_even_depth(x, rngseed=1, sample.size= 1, replace=T)
}
) 

collec

### Sample data 
sdata <- sample_data(PS)
write.csv(sdata, "/SAN/Victors_playground/Metabarcoding/sample_data.csv")
sdata <- read.csv("/SAN/Victors_playground/Metabarcoding/sample_data.csv")
names(sdata)
rpers <- read.csv("/SAN/Victors_playground/Metabarcoding/Total_reads_per_sample.csv")

sdata <- sdata[c(1,5:8)]
colnames(sdata) <- c("Position", "DNA_conc", "X260_280", "X260_230", "Species")

sdata$Position <- gsub(pattern = " ", replacement = "", x = sdata$Position)
rpers$Position <- gsub(pattern = " ", replacement = "", x = rpers$Position)

rpers <- merge(sdata, rpers, by= "Position")

##Plot DNA characteristics and total reads

library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(plm)

a1 <- ggplot(rpers, aes(x= DNA_conc, y= Total_reads,  color=Species)) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8) + 
  labs(x= "DNA concentration ng/µL", y = "Total sequencing reads") +
  scale_y_continuous(labels = scientific)+
  geom_vline(xintercept=30, linetype="dashed", color = "black")+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  #scale_x_log10()+
  #annotation_logticks()+
  geom_smooth(method = "lm", se = FALSE, col = "lightgrey") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()

rperDNA <- lm(rpers$Total_reads~rpers$DNA_conc)
summary(rperDNA)

b1 <- ggplot(rpers, aes(x= X260_280, y= Total_reads,  color=Species)) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8) + 
  geom_vline(xintercept=1.7, linetype="dashed", color = "black")+ 
  labs(x= "Purity 260/280", y = "Total sequencing reads") +
  scale_y_continuous(labels = scientific)+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  #scale_x_log10()+
  #annotation_logticks()+
  geom_smooth(method = "lm", se = FALSE, col = "lightgrey") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()

rper280 <- lm(rpers$Total_reads~rpers$X260_280)
summary(rper280)

c1 <- ggplot(rpers, aes(x= X260_230, y= Total_reads,  color=Species)) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8) + 
  geom_vline(xintercept=1.5, linetype="dashed", color = "black")+ 
  labs(x= "Purity 260/280", y = "Total sequencing reads") +
  scale_y_continuous(labels = scientific)+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  #scale_x_log10()+
  #annotation_logticks()+
  geom_smooth(method = "lm", se = FALSE, col = "lightgrey") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()

rper230 <- lm(rpers$Total_reads~rpers$X260_230)
summary(rper230)

library(gridExtra)
library(grid)

grid.arrange(a1, b1, c1, nrow= 3)

ggplot(rpers, aes(x= X260_280, y= X260_230,  color= Species)) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) + geom_hline(yintercept=1.5, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="Purity 260/280", y = "Purity 260/230") +
  theme_classic()



#####Richness and diversity 

names(STNC)[keep]<-gsub(pattern = "-", replacement = "_", x= names(STNC[keep]))
names(STNC)[keep]<-gsub(pattern = " ", replacement = "", x= names(STNC[keep]))
names(PS.l)<-names(STNC[keep])

n16S <- subset(primerInput, primerInput$Gen== "16S")
n18S <- subset(primerInput, primerInput$Gen== "18S")
n28S <- subset(primerInput, primerInput$Gen== "28S")
COI <-  subset(primerInput, primerInput$Gen== "COI")
plants <- subset(primerInput, primerInput$Target== "Viridiplantae")

setdiff(names(PS.l), as.vector(primerInput$Primer_name))


Index<- "Evenness"
vector_of_primers<- as.vector(primerInput$Primer_name) ##Change for the different groups


finalPrime <- list(list())
for(i in vector_of_primers)
{
  if(i %in% names(PS.l))
  {
    id <- which(i == names(PS.l))
    finalPrime[[i]][["richness_index"]] <- data.frame(estimate_richness(PS.l[[id]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
    Div <- finalPrime[[i]][["richness_index"]]$Shannon
    Rich <- log(finalPrime[[i]][["richness_index"]]$Observed + 1)
    evenness = Div/Rich
    finalPrime[[i]][["richness_index"]]$Evenness = evenness
    finalPrime[[i]][["number_of_samples"]] <- nrow(sample_data(PS.l[[id]]))
    mean_ind <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = mean)
    std_dev <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = sd)
    ci <- qnorm(0.95)*std_dev/sqrt(finalPrime[[i]][["number_of_samples"]])
    upper <- mean_ind+ci
    lower <- mean_ind-ci
    finalPrime[[i]][["central_tendency"]] <- data.frame(Index_name = c("Observed","Chao1", "se.chao1", "Shannon", "Simpson", "InvSimpson", "Evenness"), 
                                                        Mean = mean_ind, 
                                                        SD = std_dev,
                                                        Upper = upper, 
                                                        Lower = lower)
    finalPrime[[i]][["primer_info"]] <- data.frame(primerInput[ primerInput$Primer_name %in% i, c(6:ncol(primerInput))])
  }
}
finalPrime[[1]] <- NULL
# plots
plot_df <- data.frame()
x <- data.frame()
for(j in names(finalPrime))
{
  x <- finalPrime[[j]][["richness_index"]]
  x[1:nrow(x),8] <- j
  x[1:nrow(x),9] <- finalPrime[[j]][["primer_info"]]$Gen
  x[1:nrow(x),10] <- finalPrime[[j]][["primer_info"]]$Target
  plot_df <- rbind(plot_df, x)
}
colnames(plot_df)[8] <- c("primer_name")
colnames(plot_df)[9] <- "gene"
colnames(plot_df)[10] <- "target"

ggplot(plot_df, aes(x= primer_name, y = plot_df[,Index], fill = gene))+
  geom_boxplot(outlier.colour="black", alpha = 0.5)+
  geom_jitter(shape=1)+
  xlab("Primer name")+
  ylab(paste0(Index, " (Shannon Index/log(Observed richness +1))", collapse = ''))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1))

#pdf("AA_Primer_evaluation/Phyloseq_Test/Shannon_all.pdf", width=15, height=15)

# upsetr
require(UpSetR)
# We need the list of sample names for a given vector of primers
upset_list <- list()
for(k in vector_of_primers)
{
  upset_list[[k]] <- rownames(finalPrime[[k]][["richness_index"]])
}
upset(fromList(upset_list), order.by = "freq", mainbar.y.label = "Number of samples", sets.x.label = "Sample size", nintersects = NA, nsets = 55)



###Central tendency also worked
#stats_df <- group_by(plot_df,primer_name) %>%
#  dplyr::summarise(
#    count = n(),
#    meanShan = mean(Shannon, na.rm = TRUE),
#    sdShan = sd(Shannon, na.rm = TRUE),
#    seShan = sdShan/ sqrt(count),  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#    ciMult = qt(.95/2 + .5, n()-1),
#    ciShan = seShan * ciMult,
#    ciminShan = meanShan - ciShan,
#    cimaxShan = meanShan + ciShan,
    ## and for Simpson Index
#    meanSimp = mean(Simpson, na.rm = TRUE),
#    sdSimp = sd(Simpson, na.rm = TRUE),
#    seSimp = sdSimp / sqrt(count),  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#    ciMult = qt(.95/2 + .5, n()-1),
#    ciSimp = seSimp * ciMult,
#    ciminSimp = meanSimp - ciSimp,
#    cimaxSimp = meanSimp + ciSimp
#  )

##one way ANOVA 

#n16S <- subset(plot_df, plot_df$gene== "16S")
#n18S <- subset(plot_df, plot_df$gene== "18S")
#n28S <- subset(plot_df, plot_df$gene== "28S")
#COI <- subset(plot_df, plot_df$gene== "COI")
#plants <- subset(plot_df, plot_df$target== "Viridiplantae")
#animal <- subset(plot_df, plot_df$target== "Metazoa")


#aov(Shannon ~ primer_name, data = COI)
#summary(aov(Shannon ~ primer_name, data = COI))

#TukeyHSD(aov(aov(Shannon ~ primer_name, data = COI)))



