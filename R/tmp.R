library(phyloseq)


recomputeBioinfo <- FALSE

if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
    }
}


traits <- readRDS("input_data/traits_grouped.rds")


PSHelm <- phyloseq::subset_taxa(PS, phylum%in%c("Nematoda", "Platyhelminthes"))
PSHelmG <- phyloseq::tax_glom(PSHelm, "genus")


table(tax_table(PSHelmG)[, "genus"] %in% traits$t.genus)


## Some missing from Caro's table
tax_table(PSHelmG)[, "genus"][!tax_table(PSHelmG)[, "genus"]%in% traits$t.genus]

## Citellinema -- A trichostongyloid nematode parasite of squirrels?!

## Protostrongylus -- small ruminants, with snail i.m. host

## Neogogatea --- a Diplostomid, but otherwise NA. Nothing really to
## find about this genus

## Plagiostomum --- not a parasite! Turbellaria is the only free
## living order in Plathelminthes


write.csv(traits, "input_data/helminth_traits_.csv", row.names=FALSE)


read.csv("input_data/helminth_traits_.csv")
