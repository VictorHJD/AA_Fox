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





### IN PROGRESS ###
#Rangstufen anzeigen lassen
PSplant <- tax_glom(PSplant) #pflants zusammenfassen
tax_table(PSplant)

PSplantfamily <- tax_glom(PSplant, taxrank="family")
PSplantfamily
tax_table(PSplantfamily)


# Sequenzen aus Namen abschneiden (???) #

head(tax_table(PSFoxDiet))
tax_table(PSFoxDiet)[,"genus"]
nameOtuByTax <- function(PSFoxDiet, taxon="genus"){
  otable <- otu_table(PSFoxDiet)
  colnames(otable) <- make.unique(as.character(tax_table(PSFoxDiet)[, "genus"]))
  otable
}


## FINAL SUBSET ##

get_taxa_unique(PSFoxHigh, "phylum")
PSFoxDiet <- subset_taxa(PSFoxHigh, phylum%in%c("Arthropoda",
                                                "Mollusca", "Annelida","Streptophyta",
                                                "Chordata","Ascomycota","Basidiomycota"))

grep(pattern = "Basidiomycota", tax_table(PSFoxDiet)[,2])

tax_table(PSFoxDiet)



## SEX / AGE / CONDITION ##

table(sample_data(PSFoxDiet)$sex..m.f., sample_data(PSFoxDiet)$age)

sex <- table(sample_data(PSFoxDiet)$sex..m.f.) # turn data about sex into an object
barplot(sex, main="Distribution of sex")

table(sample_data(PSFoxDiet)$condition)

condition <- table(sample_data(PSFoxDiet)$condition)
barplot(condition, main="Distribution of condition")



## URBANISATION CATEGORIES ##

# add new column with category of urbanisation
imperviousness1000 <- table(sample_data(PSFoxDiet)$imperviousness1000)
sample_data(PSFoxDiet)$urb_factor <- ifelse(sample_data(PSFoxDiet)$imperviousness1000<25,"rural",
                                            ifelse(sample_data(PSFoxDiet)$imperviousness1000>50,"urban",
                                                   "suburban"))

# View(madeleineFullData)
length(which(sample_data(PSFoxDiet)$urb_factor=="rural")) # gives you the number of how many samples are "rural"
length(which(sample_data(PSFoxDiet)$urb_factor=="suburban"))
length(which(sample_data(PSFoxDiet)$urb_factor=="urban"))

counts <- table(sample_data(PSFoxDiet)$urb_factor)c
barplot(counts, main="Urbanisation categories")



## TIME DISTRIBUTION ##

dfforhistogram <- data.frame(data_found = sample_data(PSFoxDiet)$date_found,
                             age = sample_data(PSFoxDiet)$age_.j.a.)

dfforhistogram$monthsoftheyear <- lubridate::month(as.POSIXlt(
  dfforhistogram$data_found, format="%d/%m/%Y"), label = TRUE)

ggplot(data = dfforhistogram, aes(x = monthsoftheyear, fill = age)) +
  geom_bar(position = "dodge") +
  theme_bw() +
  scale_y_continuous(breaks=0:16)



# MAP #

library(ggplot2)
library(ggmap)

sample_data(PSFoxDiet)$location_lat <- as.numeric(as.character(sample_data(PSFoxDiet)$location_lat))
sample_data(PSFoxDiet)$location_long <- as.numeric(as.character(sample_data(PSFoxDiet)$location_long))

is.numeric(sample_data(PSFoxDiet)$location_long)

range(na.omit(sample_data(PSFoxDiet)$location_lat))
range(na.omit(sample_data(PSFoxDiet)$location_long))

area <- get_map(location = c(12, 51.5, 15, 53.5),
                source = "stamen", 
                maptype = "toner-lite", 
                zoom = 8)

map <- ggmap(area)
map

# add layers of points
map +
  theme_bw() + 
  geom_point(data = sample_data(PSFoxDiet),
             aes(x = location_long, 
                 y = location_lat))

# Pretty pretty prettier
map +
  theme_bw() + 
  geom_point(data = sample_data(PSFoxDiet),
             aes(x = location_long, 
                 y = location_lat, 
                 fill = imperviousness1000), # here put the variable you want to see
             pch = 21,
             size = 4,
             alpha = 0.8)



## ABUNDANCE BAR PLOT

# add new column with category of urbanisation
imperviousness1000 <- table(data$imperviousness1000)
data$urb_factor <- ifelse(data$imperviousness1000<25,"rural",
                          ifelse(data$imperviousness1000>50,"urban",
                                 "suburban"))

# exclude canidae
PSFox <- subset_taxa(PSFoxDiet, !family%in%"Canidae")

library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())
library(phyloseq)


#Transform into relative abundances + filter to a mean threshold
physeq2 = filter_taxa(PSFox, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3

#Turn all OTUs into class counts
glom <- tax_glom(physeq3, taxrank = 'class')
glom # should list # taxa as # class
data <- psmelt(PSFox) # create dataframe from phyloseq object
data$class <- as.character(data$class) #convert to character

#group dataframe by class, calculate median rel. abundance
library(plyr)
medians <- ddply(data, ~class, function(x) c(median=median(x$Abundance)))

#plot with condensed class
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=class))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
  facet_grid(.~urb_factor)

#compare to plot of same data, without condensed phyla
p <- ggplot(data=stom, aes(x=Sample, y=Abundance, fill=class))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))+
  facet_grid(.~urb_factor)


## ABUNDANCE BAR PLOT WITH FAMILY < 1% ABUND. IN EXTRA GROUP ##

library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())
library(phyloseq)


#Transform into relative abundances + filter to a mean threshold
physeq2 = filter_taxa(PSFoxDiet, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3

# Condense low abundance taxa into an "Other" category for the barplot
#Turn all OTUs into family counts
glom <- tax_glom(physeq3, taxrank = 'family')
glom # should list # taxa as # family
data <- psmelt(glom) # create dataframe from phyloseq object
data$family <- as.character(data$family) #convert to character

#simple way to rename family with < 1% abundance
data$family[data$Abundance < 0.01] <- "< 1% abund."

#group dataframe by family, calculate median rel. abundance
library(plyr)
medians <- ddply(data, ~family, function(x) c(median=median(x$Abundance)))

#find family whose rel. abund. is less than 1%
remainder <- medians[medians$median <= 0.01,]$family

#list of low abundance family
remainder

#change their name to "family < 1% abund."
data[data$family %in% remainder,]$family <- "family < 1% abund."
#rename family with < 1% relative abundance
data$family[data$Abundance < 0.01] <- "family < 1% abund."

#plot with condensed family into "< 1% abund" category
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=family))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
  facet_grid(.~urb_factor)

#compare to plot of same data, without condensed phyla
p <- ggplot(data=stom, aes(x=Sample, y=Abundance, fill=phylum))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))+
  facet_grid(.~urb_factor)



## JUST TOP 10 FAMILIES ##

#agglomerate OTUs at a given taxonomic level and select the 10 most abundant taxa
PSFoxfam <- tax_glom(PSFoxDiet, "family")
Fox.fam <- names(sort(taxa_sums(PSFoxfam), decreasing = TRUE)[1:10])
fam10 <- prune_taxa(Fox.fam, PSFoxDiet)


#Transform into relative abundances + filter to a mean threshold
physeq2 = filter_taxa(fam10, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3

#Turn all OTUs into family counts
glom <- tax_glom(physeq3, taxrank = 'family')
glom # should list # taxa as # family
data <- psmelt(glom) # create dataframe from phyloseq object
data$family <- as.character(data$family) #convert to character

#group dataframe by family, calculate median rel. abundance
medians <- ddply(data, ~family, function(x) c(median=median(x$Abundance)))

#plot with condensed family
library(randomcoloR)
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=family))
p + geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill=guide_legend(nrow=5))+
  facet_grid(.~urb_factor)


#export excel table
library(openxlsx)
write.xlsx(data, "/home/madeleine/foxdiet.xlsx")
write.csv(data, "/home/madeleine/foxdiet.csv")

# show all families
data$order <- as.factor(data$order)
str(data)

