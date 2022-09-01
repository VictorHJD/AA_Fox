library(SpiecEasi)
library(igraph)
library(phyloseq)

recomputeBioinfo <- FALSE

if(!exists("PS")){
    if(recomputeBioinfo){
        source("R/1_Fox_general_MA.R")
    } else {
        PS <- readRDS(file="intermediate_data/PhyloSeqCombi.Rds")
    }
}

## collapse to genus level
PSG <- phyloseq::tax_glom(PS, "genus")

PSG.target <- subset_taxa(PSG, phylum %in% c("Nematoda", "Platyhelminthes",
                                             "Chordata", "Annelida", "Arthropoda","Mollusca"))


PSG.target.bac <- subset_taxa(PSG, phylum %in% c("Nematoda", "Platyhelminthes",
                                                 "Chordata", "Annelida", "Arthropoda",
                                                 "Mollusca",
                                                 "Actinobacteria", "Bacteroidetes",
                                                 "Deferribacteres", "Firmicutes",
                                                 "Fusobacteria", "Proteobacteria",
                                                 "Spirochaetes", "Tenericutes"))



## for comparison:
naive.cor <- cor(otu_table(PSG.target), method="pearson")
fivenum(abs(naive.cor[upper.tri(naive.cor)]))

taxa_names(PSG.target) <- make.unique(tax_table(PSG.target)[, "genus"])

## taxa_names(PSHelm) <- ifelse(tax_table(PSHelm)[, "phylum"] %in%"Chordata",
##                                       as.vector(tax_table(PSHelm)[, "family", drop=T]),
##                                ifelse(tax_table(PSHelm)[, "phylum"] %in%
##                                       c("Nematoda", "Platyhelminthes"),
##                                       as.vector(tax_table(PSHelm)[, "genus", drop=T]),
##                                       as.vector(tax_table(PSHelm)[, "class", drop=T])))

## taxa_names(PSDietCollapsed) <- ifelse(tax_table(PSDietCollapsed)[, "phylum"] %in%"Chordata",
##                                       as.vector(tax_table(PSDietCollapsed)[, "family", drop=T]),
##                                ifelse(tax_table(PSDietCollapsed)[, "phylum"] %in%
##                                       c("Nematoda", "Platyhelminthes"),
##                                       as.vector(tax_table(PSDietCollapsed)[, "genus", drop=T]),
##                                       as.vector(tax_table(PSDietCollapsed)[, "class", drop=T])))


reCreateNetwork <- FALSE

if(!exists("seHelmDiet")){
    if(!reCreateNetwork){         
        seHelmDiet <- readRDS("intermediate_data/SE_networkd.Rds")
    } else {
        pargs <- list(rep.num=1000, seed=10010, ncores=20)
        seHelmDiet  <- spiec.easi(PSG.target, method="mb", pulsar.params=pargs)
    }
}

HelmDiet.mb <- adj2igraph(getRefit(seHelmDiet),
                          vertex.attr=list(name=taxa_names(PSG.target),
                                           nReads=colSums(otu_table(PSG.target))))

tax_table(PSG.target)[, "species"] <- paste0(tax_table(PSG.target)[, "genus"],
                                             " R:",  colSums(otu_table(PSG.target)),
                                             " S:",  colSums(otu_table(PSG.target)>0))

png(filename = "figures/Diet_Helm_Network.png",
        width =15, height = 15, units = "in", res= 300)
phyloseq::plot_network(HelmDiet.mb, PSG.target, type="taxa", label="species", color="phylum")
dev.off()


## Trichostrongylus is there because of a rabbit in the same sample!!!
### Trichostrongylus -- Oryctolagus



## But we also have to exclude those many taxa only in one sample
## (better even allow only those in three)

PSG.targetF <- prune_taxa(colSums(otu_table(PSG.target)>0)>3, PSG.target)

pargs <- list(rep.num=1000, seed=10010, ncores=20)
seHelmDietF  <- spiec.easi(PSG.targetF, method="mb", pulsar.params=pargs)

HelmDietF.mb <- adj2igraph(getRefit(seHelmDietF),
                          vertex.attr=list(name=taxa_names(PSG.targetF),
                                           nReads=colSums(otu_table(PSG.targetF))))

png(filename = "figures/Diet_Helm_NetworkF.png",
        width =15, height = 15, units = "in", res= 300)
phyloseq::plot_network(HelmDietF.mb, PSG.targetF, type="taxa", label="species", color="phylum")
dev.off()


######### Now with bacteria!!!! ######


PSG.targetB <- prune_taxa(colSums(otu_table(PSG.target.bac)>0)>3, PSG.target.bac)

pargs <- list(rep.num=1000, seed=10010, ncores=20)
seHelmDietB  <- spiec.easi(PSG.targetB, method="mb", pulsar.params=pargs)

HelmDietB.mb <- adj2igraph(getRefit(seHelmDietB),
                           vertex.attr=list(name=taxa_names(PSG.targetB),
                                            nReads=colSums(otu_table(PSG.targetB))))

png(filename = "figures/Bacteria_Helm_NetworkF.png",
        width =15, height = 15, units = "in", res= 300)
phyloseq::plot_network(HelmDietB.mb, PSG.targetB, type="taxa", label="species", color="phylum")
dev.off()








## Only relevant for (MA) method development!
## LET'S do this for individual Amplicons!!! 

## PS.targetF <- prune_taxa(colSums(otu_table(PS.target)>0)>3, PS.target)

## tax_table(PS.targetF)[, "species"] <- paste0(tax_table(PS.targetF)[, "genus"],
##                                              " R:",  colSums(otu_table(PS.targetF)),
##                                              " S:",  colSums(otu_table(PS.targetF)>0))

## pargs <- list(rep.num=1000, seed=10010, ncores=20)
## seHelmDietSINGF  <- spiec.easi(PS.targetF, method="mb", pulsar.params=pargs)

## HelmDietSINGF.mb <- adj2igraph(getRefit(seHelmDietSINGF),
##                                vertex.attr=list(name=taxa_names(PS.targetF),
##                                                 nReads=colSums(otu_table(PS.targetF))))

## png(filename = "figures/Diet_Helm_NetworkSINGF.png",
##         width =15, height = 15, units = "in", res= 300)
## phyloseq::plot_network(HelmDietSINGF.mb, PS.targetF, type="taxa", label="species", color="phylum")
## dev.off()


## pdf("figures/Diet_Helm_NetworkSINGF.pdf", width =15, height = 15)
## phyloseq::plot_network(HelmDietSINGF.mb, PS.targetF, type="taxa", label="species", color="phylum")
## dev.off()


## wtc <- cluster_walktrap(HelmDietSINGF.mb)

## modularity(HelmDietSINGF.mb, membership(wtc))


## pdf("figures/Diet_Helm_MODUL_net.pdf", width = 15, height = 15)
## plot(wtc,
##      HelmDietSINGF.mb,
##      label="species",
##      vertex.size=1
##      )
## dev.off()


## ### THIS IS AMAZING!!! WILL REVOLUTIONIZE MULTIAMPLICON (amplicon merging)!!!

## ### To have a better look: remove vulpes and other fuzzy things!

## ## how does diet influence Helminth occurence?
## PS.targetX <- subset_taxa(PS.targetF, !genus%in%c("Vulpes", "Procyon", "Canis"))


## tax_table(PS.targetX)[, "species"] <- paste0(tax_table(PS.targetX)[, "genus"],
##                                              " R:",  colSums(otu_table(PS.targetX)),
##                                              " S:",  colSums(otu_table(PS.targetX)>0))

## pargs <- list(rep.num=1000, seed=10010, ncores=20)
## seHelmDietSINGX  <- spiec.easi(PS.targetX, method="mb", pulsar.params=pargs)

## HelmDietSINGX.mb <-
##     adj2igraph(getRefit(seHelmDietSINGX),
##                vertex.attr=list(name=tax_table(PS.targetX)[, "species"],
##                                 nReads=colSums(otu_table(PS.targetX)),
##                                 genus=tax_table(PS.targetX)[, "genus"],
##                                 family=tax_table(PS.targetX)[, "family"],
##                                 order=tax_table(PS.targetX)[, "order"],
##                                 class=tax_table(PS.targetX)[, "class"],
##                                 phylum=tax_table(PS.targetX)[, "phylum"],
##                                 What=ifelse(tax_table(PS.targetX)[, "phylum"] %in%
##                                              c("Nematoda", "Platyhelminthes"),
##                                              "Helm", "Diet")))

## HelmDietSINGX.mb$wt.1




## png(filename = "figures/Diet_Helm_NetworkSINGX.png",
##         width =15, height = 15, units = "in", res= 300)
## phyloseq::plot_network(HelmDietSINGX.mb, PS.targetX, type="taxa", label="species", color="phylum")
## dev.off()


## pdf("figures/Diet_Helm_NetworkSINGX.pdf", width =15, height = 15)
## phyloseq::plot_network(HelmDietSINGX.mb, PS.targetX, type="taxa", label="species", color="phylum")
## dev.off()


## wtc <- cluster_walktrap(HelmDietSINGX.mb)

## modularity(HelmDietSINGX.mb, membership(wtc))


## pdf("figures/Diet_Helm_MODUL_netX.pdf", width = 15, height = 15)
## plot(wtc,
##      HelmDietSINGX.mb,
##      label="species",
##      vertex.size=1
##      )
## dev.off()

## wtc
