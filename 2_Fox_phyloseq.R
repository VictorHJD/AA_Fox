PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/PhyloSeqCombi.Rds")

nameOtuByTax <- function(ps, taxon="genus"){
    otable <- otu_table(ps)
    rownames(otable) <- make.unique(as.character(tax_table(ps)[, "genus"]))
    otable
}

head(sample_data(PS))

PSFox <- subset_samples(PS, species%in%"red fox")
PS.lFox <- lapply(PS.l, function (x) subset_samples(x, species%in%"red fox"))

ordi <- ordinate(PSFox, method="MDS", distance="bray")

plot_ordination(PSFox, ordi, color="imperviousness1000")

sample_data(PSFox)$readsum <-  sample_sums(PSFox)

plot_ordination(PSFox, ordi, color="readsum", label="readsum")


PSFoxHigh <- prune_samples(sample_sums(PSFox)>6650, PSFox)

ordi2 <- ordinate(PSFoxHigh, method="MDS", distance="bray")

plot_ordination(PSFoxHigh, ordi2, color="imperviousness1000", label="readsum")

plot_ordination(PSFoxHigh, ordi2, color="imperviousness1000", label="date_found")

ordi3 <- ordinate(subset_taxa(PSFoxHigh, !phylum%in%"Chordata"),
                  method="MDS", distance="bray")

plot_ordination(subset_taxa(PSFoxHigh, !phylum%in%"Chordata"),
                ordi3, color="imperviousness1000")

plot_ordination(subset_taxa(PSFoxHigh, !phylum%in%"Chordata"),
                ordi3, color="readsum",
                label="date_found")


library(ggplot2)

plot_richness(PSFoxHigh, "imperviousness1000", measures="Chao1") + stat_smooth()

plot_richness(PSFoxHigh, "c..ng.µl._original", measures="Chao1") + stat_smooth()

plot_richness(PSFoxHigh, "treecover1000", measures="Chao1") + stat_smooth()


## now focussing on plants
PSplant <- subset_taxa(PSFoxHigh, phylum%in%"Streptophyta")


plot_richness(PSplant, "imperviousness1000", measures="Chao1") + stat_smooth()

plot_richness(PSplant, "imperviousness1000", measures="Observed") + stat_smooth()

plot_richness(PSplant, "c..ng.µl._original", measures="Chao1") + stat_smooth()

plot_richness(PSplant, "treecover1000", measures="Chao1") + stat_smooth()


## and on invertrebrates

PSinverts <- subset_taxa(PSFoxHigh, phylum%in%c("Arthropoda",
                                                "Mollusca", "Annelida"))


plot_richness(PSinverts, "imperviousness1000", measures="Chao1") + stat_smooth()

plot_richness(PSinverts, "imperviousness1000", measures="Observed") + stat_smooth()

plot_richness(PSinverts, "c..ng.µl._original", measures="Chao1") + stat_smooth()

plot_richness(PSinverts, "treecover1000", measures="Chao1") + stat_smooth()
