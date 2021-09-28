
# read in and prepare data
Phylo <- readRDS("PhyloSeqCombi.Rds") # metabarcoding data
obj <- parse_phyloseq(Phylo)
samp_phylo <- obj$data$sample_data # just use sample data

geo <- readr::read_csv("geoCOMPLETE.csv") # sample locations
colnames(geo)[1] <- "IZW_ID" # rename column for later merging
geo$location_lat <- gsub("[^0-9\\.]", "", geo$location_lat) # covert to correct format

combi_geo <- merge(samp_phylo, geo, by = "IZW_ID", all = T) # merge sample data with sample location data
predata <- combi_geo %>% 
  select(IZW_ID,  date_found, location_lat.y, location_long.y, area, sex, age, weight_kg, `LLBB ID`, readsF) # select just relevant columns
predataNAr <- predata[complete.cases(predata[ , 10]),] # remove all rows containing mainly NAs based on column 'readsF'
 
Sys.setlocale(category = "LC_ALL", locale = "us")
imp <- readRDS("samples_sf_imp_ext_500-1000.rds") # imperviousness data
imp$date_found.y <- as_date(imp$date_found.y, format = "%d.%m.%Y", tz = "UTC") # date into date format
imp$month <- month(imp$date_found.y, label = T) # add column with month
combi_imp<- merge(predataNAr, imp, by = "IZW_ID") # merge sample data with imperviousness
prefinal <- combi_imp %>% 
  select(IZW_ID,  date_found, location_lat.y.x, location_long.y.x, area, sex, age, weight_kg, `LLBB ID`, month, imperviousness_20m, imp_500m, imp_1000m) %>% rename(lat = location_lat.y.x, long = location_long.y.x, imp_20m = imperviousness_20m ) # select and rename relevant columns


write_rds(prefinal, "./prefinal_sampledata.rds")



```


```{r}

## checking distribution of the data
prefinal %>% group_by(sex) %>% summarise(n = n()) # 75 females, 125 males
prefinal %>% group_by(age) %>% summarise(n = n()) # 183 adults, 17 juvenils
prefinal %>% group_by(area) %>% summarise(n = n()) # 147 Berlin, 53 Brandenburg
prefinal %>% group_by(month) %>% summarise(n = n()) # lack of data for summer month, especially May, June, July and August

ggplot(prefinal) +
   stat_count(aes(x = age), geom = "col") # visual illustration, just change x within aestethics
ggplot(prefinal) +
    geom_histogram(aes(x = imp_1000m)) # distribution of imperviousness values -> 1000m seems to be a good/meaningful explanatory variable




```


```{r}

## prepare parasite metabarcoding data

tax <- obj$data$tax_data
tax_sub <- tax %>% filter(phylum %in% c("Nematoda", "Platyhelminthes"))
(parasite_orders <- tax %>% filter(phylum %in% c("Nematoda", "Platyhelminthes")) %>% group_by(phylum, order) %>%  summarise(n = n()))


otu <- obj$data$otu_table
#otu_sub <- otu %>% filter(phylum %in% c("Nematoda", "Platyhelminthes"))
#View(otu)

combi_otu <- merge(otu, tax_sub, by = "otu_id", all.y = T) # merge taxa with abundance table
combi_otu_binary <- combi_otu %>%  mutate_at((3:219), funs(ifelse(. > 0, 1, .))) # change all reads into binary values, but just sample data

combi_otu_grouped <- combi_otu_binary %>% group_by(taxon_id.y) %>% summarise_at(vars(3:219), funs(sum)) # group by same otu and same values columnwise 
combi_otu_binary_grouped <- combi_otu_grouped %>%  mutate_at((2:217), funs(ifelse(. > 0, 1, .))) # change all reads into binary values again, but just sample data
taxa <- combi_otu_binary %>% group_by(taxon_id.y) %>% select(220:226) %>% distinct() # extract distinct taxa names
hopefully_final <- merge(combi_otu_binary_grouped,taxa, by = "taxon_id.y", all.x = T) # merge taxa back to binary table

write_rds(hopefully_final, "./predata_binary_MB.rds")

## adding rowsum and columnsum for prevalence and richness
hopefully_final <- hopefully_final %>% mutate(prevalence = select(., 2:218) %>% rowSums(na.rm = TRUE)) # sum of each row and therefore otu = prevalence within population
final <- hopefully_final %>% bind_rows(summarise_all(., funs(if(is.numeric(.)) sum(.) else "Total"))) # sum for each sample, total = richness per individual

readr::write_rds(final, "./final_MB_withoutsampledata.rds") # error message currently
saveRDS(final, "./final_MB_withoutsampledata.rds")

```



```{r}
## prepare metabarcoding diet data

tax <- obj$data$tax_data

## prepare diet data relevant for parasite transmission
taxa_diet <- tax %>% filter(superkingdom %in% c("Eukaryota")) %>% filter(phylum %in% c("Arthropoda", "Chordata", "Annelida", "Mollusca")) %>% filter(!order %in% c("Carnivora", "Cetacea", "Primates", "Enterogona", "Chiroptera", "Diprotodontia", "Phyllodocida")) %>% tidyr::drop_na("order") %>% mutate(diet = case_when(
    phylum == "Arthropoda" ~ "Arthropoda", 
    phylum == "Annelida" ~ "Annelida",
    phylum == "Mollusca" ~ "Mollusca",
    TRUE ~ order)) %>% mutate(diet = case_when(diet == "Siluriformes" | diet == "Syngnathiformes" | diet == "Cypriniformes" | diet == "Perciformes" | diet == "Characiformes"~ "fish", TRUE ~ diet)) %>% mutate(diet = case_when(diet == "Passeriformes" | diet == "Galliformes" | diet == "Anseriformes" | diet == "Columbiformes" ~ "bird", TRUE ~ diet))

diet_taxonID <- taxa_diet %>% select("taxon_id", "diet") %>%  distinct()## just already filtered diet (relevant taxa) and corresponding taxon_id for merging with sample data


otu <- obj$data$otu_table

combi_otu_diet <- merge(otu, diet_taxonID, by = "taxon_id", all.y = T) # merge taxa with abundance table
combi_otu_binary_diet <- combi_otu_diet %>%  mutate_at((3:219), funs(ifelse(. > 0, 1, .))) # change all reads into binary values, but just sample data
combi_otu_grouped_diet <- combi_otu_binary_diet %>% group_by(taxon_id) %>% summarise_at(vars(3:218), funs(sum)) # group by same otu and same values columnwise 
combi_otu_binary_grouped_diet <- combi_otu_grouped_diet %>%  mutate_at((2:217), funs(ifelse(. > 0, 1, .))) # change all reads into binary values again, but just sample data
diet_per_sample_binary <- merge(combi_otu_binary_grouped_diet, diet_taxonID, by = "taxon_id", all.x = T)

write_rds(diet_per_sample_binary, "./predata_diet_binary_MB.rds")

## group by diet item -> slim version (is diet item within sample - yes or no)
diet_grouped <- diet %>% group_by(diet) %>% summarise_at(vars(2:217), funs(sum)) %>% mutate_at((2:217), funs(ifelse(. > 0, 1, .)))

write_rds(diet_grouped, "./diet_binary_grouped.rds")
```



TOP TAXA - RACE PLOT
```{r}

## prevalence table 
prev_df <- hopefully_final %>% select(., -(2:218)) %>% mutate(proportion = hopefully_final$prevalence/217*100) %>% mutate_if(is.numeric, round, 2) %>% filter(prevalence != 0)   # create prevalence data frame without binary matrix, calculate the proportion of otu within population and just select cases with taxon info on genus level (157-13 cases)
prev_df_genus <- prev_df %>% tidyr::drop_na("genus") # now just 132 instead of 150 observations
prev_df_species <- prev_df %>% tidyr::drop_na("species") # 111 observations on species level

gl <- prev_df_genus %>% group_by(genus) %>% summarise(sum_genus = sum(prevalence), max_genus = max(prevalence), max_prop = max(proportion))  # prevalence on genus level; sometimes identification of same genus but different species -> unknown whether they occur in the same sample (have an impact on prevalence number and proportion), thus working with minimum proportion

gl_top <- gl %>% filter(max_prop > 20) %>% arrange(max_prop)
gl_top$genus <- c("Clonorchis - 24 %", "Oslerus - 33 %", "Uncinaria - 34 %", "Crenosoma - 37 %", "Mesocestoides - 42 %", "Eucoleus - 48 %", "Angiostrongylus - 76 %")


############## RACE PLOT #############################################

ggplot(gl_top, aes(x = reorder(genus, max_prop), y = max_prop,
    fill = genus)) + 
    geom_bar(width = 1, stat="identity") + 
    coord_polar(theta = "y") +
    xlab("") + ylab("") +
    ylim(c(0,100)) +
    ggtitle("Top Helminth Genera") +
    geom_text(data = gl_top, hjust = 1, size = 5,
              aes(x = genus, y = 0, label = genus)) +
    theme_minimal() +
    scale_fill_manual(values = c("#1b9668","#ab9a7a", "#5c5445", "#8ed0a3", "#E7B800", "#f0e472", "#2f493c")) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 18, face="bold")) 

```


```{r}

## How many foxes are infected?
t <- as.data.frame(t(final)) # reformat data frame, rows to columns
t2 <- t %>% select("V158") # just select richness column
t3 <- as.data.frame(t2[2:218, ]) # remove taxa information

IZW_ID <- row.names(t3) # extract rownames for data frame
t3$IZW_ID <- IZW_ID # add IDs 
colnames(t3) <- c("richness", "IZW_ID") # rename columns
length(which(t3$richness == 0)) # 15 samples without any detected helminth = 7%

## Impact of imperviousness and other explanatories on helminth richness, distance to city centre
rich_df <- merge(prefinal, t3, by = "IZW_ID") # combine richness values with sample data
# 52.5028889, 13.4041944 coordinates of geografical centre of Berlin
latlong <- rich_df %>% select(lat, long)
latlong$lat <- as.numeric(latlong$lat)
latlong <- SpatialPoints(latlong)
centre <- t(data.frame(c(52.5028889, 13.4041944)))
colnames(centre) <- c("lat", "long")
centre <- SpatialPoints(centre)
dist <- spDistsN1(latlong, centre, longlat = T) # calculating distance of each sample to city centre
rich_df$distance <- dist # in km
rich_df$weight_kg <- as.numeric(rich_df$weight_kg)

write_rds(rich_df, "./sampledata&richness.rds")

## test correlations between variables

GGally::ggcorr(rich_df, label = T)

## correlation between distance to city center and imperviousness?
library("ggpubr")
ggpubr::ggscatter(rich_df, x = "imp_1000m", y = "distance", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "imp", ylab = "dist")
shapiro.test(rich_df$imp_1000m) # not normal distributed
shapiro.test(as.numeric(rich_df$richness)) # not normal distributed

cor.test(rich_df$imp_1000m, rich_df$distance, method=c("spearman")) # correlation coefficient is -0.74

## dependency sex~weight
ggplot(rich_df, aes(x = sex, y = as.numeric(weight_kg))) +
    geom_boxplot(notch = T) +
    stat_summary(fun.y = mean,
        geom = "point",
        size = 3,
        color = "steelblue") +
    theme_classic() # males heavier than females

## dependency age~weight
ggplot(rich_df, aes(x = age, y = as.numeric(weight_kg))) +
    geom_boxplot(notch = T) +
    stat_summary(fun.y = mean,
        geom = "point",
        size = 3,
        color = "steelblue") +
    theme_classic() # no difference in weight

## dependency month~weight
rich_df2 <- rich_df %>% tidyr::drop_na("month")
(app1 <- ggplot(rich_df2, aes(x = month, y = as.numeric(weight_kg), fill=factor(sex))) +
    geom_boxplot() +
    scale_fill_manual(values = c("#E7B800", "#089660")) +
    ylab("weight [kg]\n") +
    guides(fill=guide_legend(title="")) +
    theme_classic()) # heavier in cold month/after autumn in weight, especcially males



rich_df3 <- rich_df %>% tidyr::drop_na("weight_kg")
plot(density(rich_df3$weight_kg)) # test distribution of weight -> gaussian
model_sex <- glm(weight_kg ~ imp_1000m + sex + sex*imp_1000m , data = rich_df, family = "gaussian")
summary(model_sex) # weigth depends on sex, but urbanisation doesn't matter



```

```{r}

## alluvial plot helminth genera ~ samples (for categorial variables)

# reframe data table -> all sample IDs for each helminth 

alluvial <- final %>% filter(prevalence != 0) %>% select(.,-c("superkingdom", "phylum", "order", "family", "species", "prevalence", "taxon_id.y")) %>% tidyr::drop_na("genus") %>% group_by(genus) %>% summarize_all(sum) %>% mutate_if(is.numeric, ~1 * (. != 0)) %>% filter(genus != "Total") %>% filter_all(any_vars(. != 0)) 

alluvial2 <- final %>% filter(prevalence != 0) %>% select(.,-c("superkingdom", "phylum", "family", "genus", "species", "prevalence", "taxon_id.y")) %>% tidyr::drop_na("order") %>% group_by(order) %>% summarize_all(sum) %>% mutate_if(is.numeric, ~1 * (. != 0)) %>% filter(order != "Total") %>% filter_all(any_vars(. != 0)) 


all.melt <- reshape2::melt(alluvial2)
#allprep <- ddply(all.melt, .(variable), transform, rescale = value)

## heatmap showing distribution of helminth genera over all samples
ggplot(all.melt, aes(variable, genus, fill = value )) +
  geom_tile() 

## sort samples along imperviousness gradient
prep_imp <- all.melt %>% dplyr::rename(IZW_ID = variable)
prep_imp2 <- left_join(prep_imp, rich_df[,c("IZW_ID", "imp_1000m")]) %>% filter(value == "1") %>% 
mutate(imperviousness = case_when(
    .$imp_1000m < 25 ~ "low",
    .$imp_1000m > 50 ~ "high",
    TRUE ~ "medium")) %>% mutate(helminth = case_when(
      .$order == "Enoplida" | .$order == "Plagiorchiida" | .$order == "Plectida" | .$order == "Prolecithophora" | .$order == "Rhabdocoela" | .$order == "Rhinebothriidea"~ "other (6)",
      TRUE ~ .$order)) 
  
  
 


ggplot(prep_imp2, aes(axis1 = order, axis2 = imperviousness, y = 1)) +
  scale_x_discrete(limits = c("order", "imp_1000m")) +
  geom_alluvium(aes(fill = order), show.legend = F) +
  geom_stratum(stat = "stratum", width = 1/2) + 
  geom_text(stat = "stratum", infer.label = TRUE) +
  theme_void()


prep_imp2$imperviousness<- factor(prep_imp2$imperviousness,levels = c("low", "medium", "high"))
prep_imp2$helminth<- factor(prep_imp2$helminth,levels = c("Opisthorchiida", "Cyclophyllidea", "Rhabditida", "other (6)", "Strongylida", "Trichinellida", "Strigeidida"))

prep_imp2 %>% 
ggplot(aes(y = 1, axis1 = helminth, axis2 = imperviousness)) +
  geom_flow(aes(fill = helminth), show.legend = F) +
  scale_fill_manual(values = c("#1b9668","#f0e472", "#5c5445", "#8ed0a3", "#E7B800", "#ab9a7a", "#2f493c")) +
  geom_stratum(color = "grey30", fill = "#f6f3ef") + 
  #geom_text(stat = "stratum", infer.label = TRUE, size = 2) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), fill = "white", size = 6) +
  theme_void()



```


```{r}

## generalised linear model to test impact of urbanisation on helminth richness

library(stats)
 # glm because richness data not normal distributed
m2 <- glm(as.numeric(richness) ~ as.numeric(imp_1000m) + as.character(sex), data = rich_df)
summary(m2)

rich_mod <- rich_df %>% filter(sex == "female")
m3 <- glm(as.numeric(richness) ~ as.numeric(imp_1000m) + as.numeric(weight_kg), data = rich_mod)
summary(m3)
cor.test(rich_mod$weight_kg, as.numeric(rich_mod$richness), method=c("spearman")) # correlation coefficient is 0.06


ggplot(rich_df, aes(x = sex, y = as.numeric(richness))) +
    geom_boxplot(notch = T) +
    stat_summary(fun.y = mean,
        geom = "point",
        size = 3,
        color = "steelblue") +
    theme_classic()

ggplot(rich_mod, aes(x = weight_kg, y = as.numeric(richness))) +
  geom_point() +
  theme_classic()

ggplot(rich_df, aes(x = distance, y = as.numeric(richness))) +
  geom_point() +
  theme_classic()

ggplot(rich_df, aes(x = weight_kg, y = as.numeric(richness))) +
  geom_point() +
  theme_classic()

## nice bubble plot with multiple variables
(b <- ggplot(rich_df, aes(x = imp_1000m, y = as.numeric(richness))) +
  geom_point(aes(color = sex, size = as.numeric(weight_kg), alpha = 0.5)) +
  scale_color_manual(values = c("#E7B800", "#089660")) +
  scale_size(range = c(0.5, 12)) +
  theme_minimal() +
  labs(x = "degree of imperviousness [%]", y = "helminth richness\n") +
  theme(legend.position="bottom",
        legend.box = "vertical", 
        axis.title.x = element_text(size = 17, face = "bold"),
        axis.title.y = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 16),
        legend.title = element_text(color = "#5c564b", size = 14, face = "bold"),
        legend.text = element_text(color = "#5c564b", size = 14)) +
  guides(alpha = F) +
  guides(size = guide_legend(override.aes = list(shape = 1))) +
  scale_size(name = "weight [kg]"))
 

library("ggExtra")
ggMarginal(b, type = "density")


c <- ggplot(rich_df, aes(x = distance, y = as.numeric(richness))) +
  geom_point(aes(color = sex, size = as.numeric(weight_kg), alpha = 0.5)) +
  scale_color_manual(values = c("#E7B800", "#089660")) +
  scale_size(range = c(0.5, 12)) +
  theme_minimal() +
  labs(x = "distance to city centre [km]", y = "helminth richness") +
  theme(legend.position="bottom",
        legend.box = "vertical", 
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        legend.title = element_text(color = "#5c564b", size = 12, face = "bold"),
        legend.text = element_text(color = "#5c564b", size = 12)) +
  guides(alpha = F) +
  guides(size = guide_legend(override.aes = list(shape = 1))) +
  scale_size(name = "weight [kg]")
 

#library("ggExtra")
ggMarginal(c, type = "density")


```


#############
rarefaction and extrapolation of helminth using iNEXT, calculate Jaccard index (eveness of categories, here low, medium and high imperviousness)

```{r}

predata_binary_MB <- readRDS( "./genus_binary_grouped.rds")
sample_data <- readRDS("./sampledata&richness.rds")

### creating master table

diet2 <- diet %>% t() %>% as.data.frame(stringsAsFactors = F) %>% janitor::row_to_names(row_number = 1) %>%  tibble::rownames_to_column("IZW_ID")

###  group binary_MB data into 3 groups: low (<25%), medium (25-50%) and high (>50%) imperviousness

# t <- predata_binary_MB %>%  
#   t() %>%
#   as.data.frame(stringsAsFactors = F) %>%
#   tibble::rownames_to_column("IZW_ID") %>%
  # `colnames<-`(.[1,]) %>%
  # .[-1,] %>%
  # `rownames<-`(NULL) %>% 
  # slice(1:217,) # rows to columns and just keep presence-abscence data

# colnames(t)[1] <- "IZW_ID"

# master <- merge(sample_data, t, by = "IZW_ID")
# master2 <- merge(imp_cat, diet2, by = "IZW_ID")
# write_rds(master2, "./sampledata&OTU&diet.rds")
  
# imp <- sample_data %>% select(IZW_ID, imp_1000m)
# combi <- merge(master, imp, by = "IZW_ID", all = F)
# imp_cat <- master %>% mutate(imperv_cat = case_when(
#     imp_1000m < 25 ~ "low", 
#     imp_1000m < 50 ~ "medium",
#     imp_1000m > 50 ~ "high"
#    ))


OTU_low <- imp_cat %>% filter(imperv_cat == "low") %>% select(.,-c("imp_1000m", "imperv_cat")) %>% column_to_rownames("IZW_ID") %>% t()
class(OTU_low) <- "numeric"

OTU_mid <- imp_cat %>% filter(imperv_cat == "medium") %>% select(.,-c("imp_1000m", "imperv_cat")) %>% column_to_rownames("IZW_ID") %>% t() 
class(OTU_mid) <- "numeric"

OTU_high <- imp_cat %>% filter(imperv_cat == "high") %>% select(.,-c("imp_1000m", "imperv_cat")) %>% column_to_rownames("IZW_ID") %>% t()
class(OTU_high) <- "numeric"


OTU_all <- list(OTU_low, OTU_mid, OTU_high)
names(OTU_all) <- c("low", "mid", "high")


### diversity analysis

DataInfo(OTU_low)

# first exploration
OTU_inext <- iNEXT(OTU_all, q =0, datatype = "incidence_raw")
OTU_inext$AsyEst # these are the asyntotic estimates (the maximum that the value will reach with infinite samples)
OTU_inext$iNextEst
# make categories comparable using the same sample size
# Maximun extrapolation is double the minimum observed sample serialize(High = 54 ind. 
# Max extrapolation: 54 x 2 + 232)
OTU_est_inext <- estimateD(OTU_all, datatype = "incidence_raw", base = "size", 
                           level = 232)
OTU_est_inext[OTU_est_inext$order ==0,] # qD is the estimated value of order 0 = especies richness, qD.LCL and qD.UCL are the lower and upper confidence intervals.

ggiNEXT(OTU_inext)

```

