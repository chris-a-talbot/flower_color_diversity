source("helper_raster.R")

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

# Load required packages
p_load("data.table", "stringr", "dplyr", 
       "parallel", "sf", "terra", "tools", "rangeBuilder")

core_num = detectCores() - 1

# Get the flowers masterfile
flowers = fread("./data/better_flowers.csv")

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get the color trait data
traits_color = as.data.frame(flowers[,7:15])
rownames(traits_color) = species_list

species_ranges = get_species_ranges(get_shapefile_names(), cores=1, debug=TRUE)

get_monthly_rasters_from_ranges(species_ranges, traits_pheno, 
                                c(-92.5, -67.5, 38, 48), res=0.09, write=TRUE, 
                                cores=1, months=c("apr", "may", "jun", "jul", "aug", "sep"))
