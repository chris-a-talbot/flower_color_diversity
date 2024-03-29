### STEP 1: Get all necessary data
### Packages, shapefiles, color traits, phenology, and phylogeny

# Pre-made R functions for handling rasters & presence/absence matrices
source("helper_raster.R")
source("helper_ecdf.R")
source("helper_metrics.R")

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

p_load("data.table", "stringr", "dplyr", "tools", "parallel",
       "terra", "abdiv")

# "ape", "phytools"

core_num = detectCores() - 1

# Get the flowers masterfile
flowers = fread("./data/better_flowers.csv")

# Get the pre-selected cells to analyze
best_indices = readRDS("./data/indices.rds")

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get the color trait data
traits_color = as.data.frame(flowers[,7:15])
rownames(traits_color) = species_list

# # Get the phylogenetic data
# tree_path = "./data/flowers_tree.txt"
# tree = read.newick(tree_path)
# tree$tip.label = str_replace(word(tree$tip.label, 1, 2, sep="_"), "_", " ")
# gen_dist_mat = cophenetic.phylo(tree)
# species_names = tree$tip.label

# Get the presence/absence data
months = c("feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov")
pa_matrices_by_month = get_month_pa_matrix_from_tif(months, res=0.09)
cell_num = nrow(pa_matrices_by_month[[1]])

# Creates color trait dataframes for each month
traits_color_by_month <- lapply(names(pa_matrices_by_month), function(month) {
  pa_matrix <- pa_matrices_by_month[[month]]
  
  # Extract species names from the layer names of the stack
  month_species_present <- colnames(pa_matrix)
  
  # Filter color_traits to include only those species
  # Ensure row names are indeed the species names in color_traits
  return(traits_color[month_species_present, , drop = FALSE])
}) 
names(traits_color_by_month) = months

# # Creates a genetic distance matrix for each month
# gen_dist_by_month <- lapply(months, function(month) {
#   pa_matrix <- pa_matrices_by_month[[month]]
#   
#   # Extract species names from the layer names of the stack
#   month_species_present <- colnames(pa_matrix)
#   
#   month_tree = drop.tip(phy = tree, tip = setdiff(tree$tip.label,
#                                                   month_species_present))
#   return(cophenetic.phylo(month_tree))
# })
# 
# tree_by_month = lapply(1:length(months), function(month) {
#   month_tree = drop.tip(phy = tree, tip = setdiff(tree$tip.label,
#                                                   colnames(pa_matrices_by_month[[month]])))
#   return(month_tree)
# })

summed_traits_lists = mclapply(months, function(month) {
  month_pa_matrix = pa_matrices_by_month[[month]]
  trait_list = get_color_distributions(month_pa_matrix, indices=best_indices, 
                                       colors=traits_color)
  return(trait_list)
}, mc.cores=core_num)
names(summed_traits_lists) = months

simpsons_es = mclapply(months, function(month) {
  traits = summed_traits_lists[[month]]
  return(get_simpsons_evenness(traits))
}, mc.cores=core_num)
names(simpsons_es) = months

summed_traits_lists_null = mclapply(months, function(month) {
  lists = lapply(1:100, function(i) {
    month_pa_matrix = pa_matrices_by_month[[month]]
    trait_list = get_color_distributions(month_pa_matrix, indices=best_indices, 
                                         colors=traits_color, null_hypothesis=TRUE)
    return(trait_list)
  })
  names(lists) = best_indices
  return(lists)
}, mc.cores=core_num)
names(summed_traits_lists_null) = months

simpsons_es_null = mclapply(months, function(month) {
  traits = summed_traits_lists_null[[month]]
  es_month = lapply(1:length(traits), function(i) {
    es_unit = lapply(1:100, function(j) {
      return(get_simpsons_evenness(traits[[i]][j,]))
    })
  })
  names(es_month) = best_indices
  return(es_month)
}, mc.cores=core_num)
names(simpsons_es_null) = months

get_cell_distributions = function(null_distributions) {
  cell_distributions = list()
  
  for(i in 1:length(null_distributions)) {
    # Temporary list to store the FunRao values for the current overarching item
    tempResults <- list()
    
    cell_num = length(null_distributions[[i]][[1]])
    
    # Initialize lists to collect FunRao values for each label
    for(j in 1:cell_num) {
      tempResults[[as.character(j)]] <- vector()
    }
    
    # Iterate over each of the 100 items within the current overarching item
    for(item in null_distributions[[i]]) {
      # Access the FunRao list

      # Iterate over each integer label in the FunRao list
      for(label in 1:length(item)) {
        # Append the value to the appropriate list in tempResults
        tempResults[[as.character(label)]] <- c(tempResults[[as.character(label)]], item[[label]])
      }
    }
    
    # Add the compiled FunRao values for the current overarching item to the results list
    cell_distributions[[as.character(i)]] <- tempResults
  }
  
  return(cell_distributions)
}

# Get the list of 100 null distributions for each site for each month
cell_distributions = mclapply(months, function(month) {
  distributions = get_cell_distributions(simpsons_es_null[[month]])
  names(distributions) = best_indices
  return(distributions)
}, mc.cores=core_num)
names(cell_distributions) = months

# Save all of this for later use
save.image("./RData/simpsons.RData")
