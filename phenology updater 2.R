### STEP 1: Get all necessary data
### Packages, shapefiles, color traits, phenology, and phylogeny

# Pre-made R functions for handling rasters & presence/absence matrices
source("helper_raster.R")
source("helper_ecdf.R")
source("helper_fu.R")
source("helper_fd.R")

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

# Load required packages
p_load("ape", "picante", "data.table", "stringr", "purrr", "dplyr", 
       "parallel", "sf", "terra", "tools")

core_num = detectCores() - 1

# Get the flowers masterfile
flowers = fread("./data/flowers_master.csv")

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get the color trait data
traits_color = as.data.frame(flowers[,7:15])
rownames(traits_color) = species_list

species_list_present = word(file_path_sans_ext(get_shapefile_names(location="./shp/pts")), 1, sep="_")

phenologies = mclapply(species_list_present, function(species) {
  phenology = fread(paste0("./data/phenology/", species, ".csv"))
  print(paste("Read phenology of", species))
  return(phenology)
}, mc.cores=core_num)

modified_phenologies <- mclapply(phenologies, function(dt) {
  # Check if the data.table is empty (0x0)
  if (nrow(dt) == 0) {
    # Replace empty data.tables with NA
    return(NA)
  } else {
    # Apply the filtering for non-empty data.tables
    output = dt[longitude >= -92.5 & longitude <= -67.5 & latitude >= 38 & latitude <= 48,]
    if(nrow(output) < 2) {
      return(NA)
    } else {
      return(output)
    }
  }
}, mc.cores=core_num )

# Correctly identifying NA 'data.table' objects or NA list elements
na_indices <- which(sapply(modified_phenologies, function(x) is.null(x) || (is.na(x) && !is.data.frame(x))))

# Removing items at these indices from both lists
modified_phenologies <- modified_phenologies[-na_indices]
species_list_present <- species_list_present[-na_indices]

flowers_present = flowers[word(taxon, 1, 2) %in% species_list_present]

flowers_modified = mclapply(1:length(species_list_present), function(species_index) {
  species_name = species_list_present[species_index]
  obs = modified_phenologies[species_index][[1]]

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

  row = flowers_present[word(taxon, 1, 2) == species_name]
  modified_row = lapply(1:12, function(month) {
    month_index = month+15
    row[,month_index] = (month %in% months)
  })
  row$inat_phen = TRUE
  row_p1 = ncol(row)+1
  row[,"ob_num"] = nrow(obs)
  return(row)
}, mc.cores=core_num)
flowers_modified = rbindlist(flowers_modified)
fwrite(flowers_modified, "./data/better_flowers.csv")