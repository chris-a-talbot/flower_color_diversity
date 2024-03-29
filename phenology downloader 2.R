library(data.table)
library(vegan)
library(dplyr)
library(gtools)
library(rinat)
library(stringr)
library(tidyverse)

# Files
wd = "G:/My Drive/Documents/Research/Weber/honors-thesis"
flowers_master_file = "data/flowers_master.csv"

# Set working directory
setwd(wd)

# Get a data.table of all flower data
flowers_master = fread(flowers_master_file)[order(species),]
num_species = nrow(flowers_master)
flowers_master = fread(flowers_master_file)[order(species),]
species_list_taxon = word(flowers_master$taxon, 1, 2)
species_list_taxon_3 = word(flowers_master$taxon, 1, 3)
species_list_name = word(flowers_master$species, 1, 2)
tree_path = "./data/flowers_tree.txt"
tree = read.newick(tree_path)
species_names = tree[["tip.label"]]
num_species = nrow(flowers_master)

last_i = 1
for(i in last_i:num_species) {
  species_name = word(flowers_master$taxon[i], 1, 2)
  if(i %in% c(50, 428, 606)) {
    species_name = substr(species_name,1,nchar(species_name)-2)
  }
  obs = get_inat_obs(taxon_name=species_name, annotation=c(12,13), maxresults=10000, bounds=c(38, -92.5, 48, -67.5), quality="research", geo=TRUE)
  if(nrow(obs) <= 3) { next }
  write.csv(obs, paste0("data/phenology/", species_name, ".csv"))
  
  dates = word(obs$datetime, 1) %>% str_split("-")
  months = as.data.table(as.numeric(sapply(dates, "[[", 2)))
  month_num = nrow(months)
  months_unique = as.list(unique(months))[[1]]
  for(j in 1:length(months_unique)) {
    val = months_unique[j]
    n = nrow(months[V1 == val])
    percent = n/month_num
    if(percent <= 0.125) { months = months[V1 != val] }
  }
  months = as.list(unique(months))[[1]]
  
  flowers_master[i, "jan"] = (1 %in% months)
  flowers_master[i, "feb"] = (2 %in% months)
  flowers_master[i, "mar"] = (3 %in% months)
  flowers_master[i, "apr"] = (4 %in% months)
  flowers_master[i, "may"] = (5 %in% months)
  flowers_master[i, "jun"] = (6 %in% months)
  flowers_master[i, "jul"] = (7 %in% months)
  flowers_master[i, "aug"] = (8 %in% months)
  flowers_master[i, "sep"] = (9 %in% months)
  flowers_master[i, "oct"] = (10 %in% months)
  flowers_master[i, "nov"] = (11 %in% months)
  flowers_master[i, "dec"] = (12 %in% months)
  flowers_master[i, "inat_phen"] = TRUE
  fwrite(flowers_master, flowers_master_file)
  pct_done = round((i/1096)*100, digits=2)
  cat(paste0(i, "/1096: ", pct_done, "% complete. \n"))
}

write.csv(flowers_master, flowers_master_file)

species_median_days = data.frame()
for(i in 1:num_species) {
  species_name = word(flowers_master$taxon[i], 1, 2)
  if(i %in% c(50, 428, 606)) {
    species_name = substr(species_name,1,nchar(species_name)-2)
  }
  obs = try(fread(paste0("data/phenology/", species_name, ".csv")))
  if(typeof(obs) != "list") next
  median_day = median(as.numeric(strftime(obs$datetime, format="%j")))
  taxon_name = str_detect(species_names, gsub('\\)', "-", gsub('\\(', "-", gsub(" ", "_", species_list_taxon_3[[i]])))) %>%
    keep(species_names, .)
  species_median_days = rbind(species_median_days, c(taxon_name, median_day))
}
colnames(species_median_days) = c("taxon", "median_day")
species_median_days = (as.data.table(species_median_days))[!is.na(median_day)]
species_median_days_final = as.data.frame(species_median_days[,-1])
rownames(species_median_days_final) = species_median_days$taxon
species_median_days_final$median_day = as.numeric(unlist(species_median_days_final$median_day))
