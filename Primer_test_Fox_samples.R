## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
library("MultiAmplicon")
## devtools::install_github("benjjneb/dada2", force= T)


###################Primer test fox#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Metabarcoding/Primer_Test_Fox" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

#samples <- gsub("_R1_001.fastq\\.gz", "\\1", basename(fastqF))

samples <- gsub("_S[0-9]*_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF)) #works but still is necessary to take out the first number
samples <- gsub("[0-9]*-Ch", "Ch", basename(samples))
samples

#Extra step in the pipeline: quality plots of the reads 
plotQualityProfile(fastqF[[1]])
plotQualityProfile(fastqF[[2]])
plotQualityProfile(fastqR[[1]])
plotQualityProfile(fastqR[[2]])

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
if(sum(file.exists(fastqF)) -
   sum(file.exists(filtFs)) > 50){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(200,200), minLen=c(200,200), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples

files <- PairedReadFileSet(filtFs, filtRs)
files
#> Warning in validityMethod(object): 
#> file ~/filtered_sra//SRR5569127_F_filt.fastq.gz does not exist on your system
#> file ~/filtered_sra//SRR5569127_R_filt.fastq.gz does not exist on your system

#Preparation of primer file 

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/primer_file_foxes.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primers <- PrimerPairsSet(primerF, primerR)

#Multi amplicon pipeline
MA <- MultiAmplicon(primers, files)

filedir <- "/SAN/Victors_playground/Metabarcoding/stratified_files"
if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
MA <- sortAmplicons(MA, n=1e+05, filedir=filedir)

#knitr::kable(getRawCounts(MA))

errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                     verbose=0)
#plotErrors(errF, nominalQ=TRUE)
errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                    verbose=0)
#plotErrors(errR, nominalQ=TRUE)

MA <- derepMulti(MA, mc.cores=4) 
MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                verbose=0)
head(MA) ##56 samples


MA <- mergeMulti(MA) ## with multiple cores MA <- mergeMulti(MA, mc.cores=20)

head(MA)
rownames(MA)

rowSums(getRawCounts(MA))
foo=rowSums(getRawCounts(MA))
foo[order(foo)]

#dput(foo, "/SAN/Victors_playground/Metabarcoding/rawcounts_per_primerpair.txt") save raw counts per primer pair 

propMerged <- MultiAmplicon::calcPropMerged(MA[c(1:6,8:21,23:75),])

head(propMerged)
(names(propMerged)) ## In this step there are 73 primer pairs 

summary(propMerged)

table(propMerged<0.8) ## For the fox test just 8 amplicons could be merged retaining over 80% of the sequence

MA <- mergeMulti(MA[c(1:6, 8:21, 23:75),], justConcatenate=propMerged<0.8) ##Works

rownames(MA)

MA <- makeSequenceTableMulti(MA) ##Bugging AGAIN, BUT primer pair 32 ("Plat-diploCOX1F-89_F.Plat-diploCOX1R-90_R") is generating the problem

MA <- makeSequenceTableMulti(MA[c(1:31,33:73)]) ##Working 

MA <- removeChimeraMulti(MA)


tracking <- getPipelineSummary(MA)

plotPipelineSummary(tracking) 
library(ggplot2)
plotPipelineSummary(tracking) + scale_y_log10()

###Extract sequences to do taxonomic assignment 

STNC <- getSequenceTableNoChime(MA)

sequences <- unlist(lapply(STNC, colnames))
names(sequences) <- paste0("asv_", 1:length(sequences))

Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
                            "/SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta")

# clusters <-
plotAmpliconNumbers(MA) 

###BLAST
#blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Metabarcoding/FoxTest_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -max_target_seqs 5 -out /SAN/Victors_playground/Metabarcoding/Alignment_fox.asn
#blast_formatter -archive /SAN/Victors_playground/Metabarcoding/Alignment_fox.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > |/SAN/Victors_playground/Metabarcoding/Alignment_fox.blttax

###Read blast result 
## we read that ouput into R blast <-
blast <- read.csv("/SAN/Victors_playground/Metabarcoding/Alignment_fox.blttax", header=FALSE)

names(blast) <- c("query", "subject", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "staxid")

library(taxonomizr)
library(taxize)

taxaNodes <- read.nodes("/SAN/db/taxonomy/nodes.dmp")
taxaNames <- read.names("/SAN/db/taxonomy/names.dmp")

blast.tax <- getTaxonomy(blast$staxid,taxaNodes,taxaNames)
blast.tax <- as.data.frame(blast.tax)
rownames(blast.tax) <- NULL

## saveRDS(blast.tax, file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
## blast.tax <- readRDS(file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")

blast <- cbind(blast, blast.tax)

library(data.table)

blt <- as.data.table(blast) 
## blt <- as.data.frame.table(blast)

head(blt)

blt <- blt[,.(bitsum=sum(bitscore),
              superkingdom, phylum, class, order, family, genus, species),
           by=c("query", "subject")]

blt <- unique(blt)
blt
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

keep <- unlist(lapply(annot.list, nrow))>0

annot.list <- annot.list[keep]
STNC <- STNC[keep]

cbind(cumsum(unlist(lapply(annot.list, nrow))), cumsum(unlist(lapply(STNC, ncol))))

## Tabulate Phyla for each amplicon
lapply(annot.list, function (x) table(x[, "phylum"]))

tabulate.genera <- function(tax, subset){
  t <- tax[phylum%in%subset, ]
  if(!is.null(ncol(t))){
    table(t[,genus])
  } else {NULL}
}

phylalist <- lapply(annot.list, function (x) table(x[, "phylum"]))
write.csv(as.data.table(phylalist), "/SAN/Victors_playground/Metabarcoding/phyla_list_fox.csv")

#phylalist <- as.data.table(phylalist)

## Tabulate by specific phylum
lapply(annot.list, function (x) tabulate.genera(x,  "Chordata"))
lapply(annot.list, function (x) tabulate.genera(x,  "Nematoda"))
lapply(annot.list, function (x) tabulate.genera(x,  "Apicomplexa"))
lapply(annot.list, function (x) tabulate.genera(x,  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.genera(x,  "Streptophyta"))


### all.annot <- Reduce(rbind, annot.list)

library(phyloseq)

p.left <- names(annot.list)
length(STNC[p.left])

annot.list <- lapply(seq_along(annot.list), function(i){
  newannot <- as.matrix(annot.list[[i]])
  taxnames <- colnames(STNC[[i]])
  rownames(newannot) <- taxnames
  newannot
})


PS.l <- lapply(seq_along(STNC), function(i){
  phyloseq(otu_table(STNC[[i]], taxa_are_rows=FALSE),
           ##  sample_data(samples.long[rownames(ALL), ]),
           tax_table(annot.list[[i]]))
})


sumSeqByTax <- function (Phy, tax) {
  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
  counts$asvCount <- as.numeric(as.character(counts$asvCount))
  tapply(counts$asvCount, counts[, tax], sum)
}

readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
names(readNumByPhylum) <- names(STNC)


readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
names(readNumByGenus) <- names(STNC)

write.csv(as.data.table(names(readNumByGenus)), "/SAN/Victors_playground/Metabarcoding/genera_list_fox.csv")

