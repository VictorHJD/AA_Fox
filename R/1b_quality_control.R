library(phyloseq)
library(ggplot2)
library(tidyverse)

## The "bioinformatics" script "1_Fox_general_MA.R" has finer controls
## within itself for now. But we can still read it's single final
## output object.
recomputeBioinfo <- FALSE

if(!exists("PSG")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PSG <- readRDS(file="intermediate_data/PhyloSeqGenus.Rds")
    }
}

### we do this without the "b" post-fixes in sample names

### The samples taken for sequencing ############################
sample.data <- readRDS("intermediate_data/Fox_data_envir.RDS")
sample.data$ID <- paste("Fox", gsub(" |b", "", sample.data$IZW_ID))


## The sampples producing good sequencing data ###################
SD <- sample_data(PSG)
SD$ID <- gsub("b", "", SD$Site)


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

## merge the Metadata and the DNA data
### 473 match                               <<- RESULT, extracted DNA!!
DNA <- merge(DNA, Meta, by="ID")

DNA$sequenced <- DNA$ID%in%sample.data$ID
table(DNA$sequenced)
## 226                                       <<- RESULT, sequenced!

DNA$analysed <- DNA$ID%in%SD$ID
table(DNA$analysed)

## only te sequenced
SEQ <- DNA[DNA$sequenced,]

#### confirming: sequenced were the better samples
tapply(DNA$ng.ul, DNA$sequenced, median)
tapply(DNA$ng.ul, DNA$sequenced, mean)


## lower DNA quality sampels tend to not be analysed futher... they
## failed in sequencing and were screened
tapply(SEQ$ng.ul, SEQ$analysed, median)
tapply(SEQ$ng.ul, SEQ$analysed, mean)

DNA$age..j.a.[DNA$age..j.a.%in%"juvenil"] <- "juvenile"


############## OVERALL EXTRACTION #############################################

summary(lm(ng.ul~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=DNA))
### lower quantity of DNA for juvenile and heavier foxes, much lower
### in Brandenburg!!!

summary(lm(X260.280~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=DNA))
### different (lower) quality in Brandenburg

summary(lm(X260.230~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=DNA))
### different (lower) quality in Brandenburg

############## WAS BIASED #####################################################




############## SEQUENCED SAMPLES  #############################################

summary(lm(ng.ul~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ))
### lower DNA for juvenile and heavier foxes, much lower in Brandenburg!!!

summary(lm(X260.280~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ))
### different (lower) quality in Brandenburg

summary(lm(X260.230~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ))

############## STILL BIASED #####################################################



############## ANALYSED SAMPLES  #############################################

summary(lm(ng.ul~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ[SEQ$analysed,]))
### PROBLEM FIXED !!!!!!!!!!!!!

summary(lm(X260.280~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ[SEQ$analysed,]))
### Well... mhmmm 

summary(lm(X260.230~I(as.numeric(weight..kg.))+condition+sex..m.f.+ age..j.a.+area,
           data=SEQ[SEQ$analysed,]))
### Well... mhmmm 

############## BIAS is REMOVED in processing ############## ##############

