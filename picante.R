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
p_load("phangorn", "phytools", "ape", "picante", "data.table", "stringr", "purrr", "dplyr", 
       "parallel", "sf", "terra", "rangeBuilder", "tools", "BiocManager")

core_num = detectCores() - 1

# Get the flowers masterfile
flowers = fread("./data/better_flowers.csv")

# Get the list of indices
indices = readRDS("./data/indices.RDS")

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get the color trait data
traits_color = as.data.frame(flowers[,7:15])
rownames(traits_color) = species_list

# Get the phylogenetic data
tree_path = "./data/flowers_tree.txt"
tree = read.newick(tree_path)
tree$tip.label = str_replace(word(tree$tip.label, 1, 2, sep="_"), "_", " ")
gen_dist_mat = cophenetic.phylo(tree)
species_names = tree$tip.label

# Get the presence/absence data
months = c("apr", "may", "jun", "jul", "aug", "sep")
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

# Creates a genetic distance matrix for each month
gen_dist_by_month <- lapply(months, function(month) {
  pa_matrix <- pa_matrices_by_month[[month]]
  
  # Extract species names from the layer names of the stack
  month_species_present <- colnames(pa_matrix)
  
  month_tree = drop.tip(phy = tree, tip = setdiff(tree$tip.label,
                                                  month_species_present))
  return(cophenetic.phylo(month_tree))
})

tree_by_month = lapply(1:length(months), function(month) {
  month_tree = drop.tip(phy = tree, tip = setdiff(tree$tip.label,
                                            colnames(pa_matrices_by_month[[month]])))
  return(month_tree)
})

# Phylogenetic dispersion tests
pd.results = mclapply(1:length(months), function(month) {
  pa_matrix_month = pa_matrices_by_month[[month]][indices,]
  tree_month = tree_by_month[[month]] %>% drop.tip(tip = setdiff(tree$tip.label, 
                                                                 colnames(pa_matrix_month)))
  pd.result = pd(pa_matrix_month, tree_month, include.root = TRUE)
  return(pd.result)
})
names(pd.results) = months

ses.pd.results = mclapply(1:length(months), function(month) {
  pa_matrix_month = pa_matrices_by_month[[month]][indices,]
  tree_month = tree_by_month[[month]] %>% drop.tip(tip = setdiff(tree$tip.label, 
                                                                 colnames(pa_matrix_month)))
  ses.pd.result = ses.pd(pa_matrix_month, tree_month, null.model="taxa.labels", runs=99)
  return(ses.pd.result)
}, mc.cores=1)
names(ses.pd.results) = months

save.image("./RData/PD.RData")

ses.mpd.results = mclapply(1:length(months), function(month) {
  pa_matrix_month = pa_matrices_by_month[[month]][indices,]
  tree_month = tree_by_month[[month]] %>% drop.tip(tip = setdiff(tree$tip.label, 
                                                                 colnames(pa_matrix_month)))
  gen_dist = cophenetic(tree_month)
  ses.mpd.result = ses.mpd(pa_matrix_month, gen_dist, null.model="taxa.labels", runs=99)
  return(ses.mpd.result)
})
names(ses.mpd.results) = months

save.image("./RData/PD.RData")

ses.mntd.results = mclapply(1:length(months), function(month) {
  pa_matrix_month = pa_matrices_by_month[[month]][indices,]
  tree_month = tree_by_month[[month]] %>% drop.tip(tip = setdiff(tree$tip.label, 
                                                                 colnames(pa_matrix_month)))
  gen_dist = cophenetic(tree_month)
  ses.mpd.result = ses.mntd(pa_matrix_month, gen_dist, null.model="taxa.labels", runs=99)
  return(ses.mpd.result)
})
names(ses.mntd.results) = months

save.image("./RData/PD.RData")











# 
# 
# 
# 
# 
# 
# # For trait tests
# colors = data.table(taxon = flowers_community$taxon, color = flowers_community$color_desc)
# color_list = c()
# for(i in 1:length(unique(flowers_community$color_desc))) {
#   colors[color == unique(flowers_community$color_desc)[i]]$color = as.numeric(i)
#   color_list = c(color_list, unique(flowers_community$color_desc)[i])
# }
# species_names = colors$taxon
# color_trait = colors[,-1]
# color_trait$color = as.numeric(unlist(color_trait$color))
# color_trait = as.data.frame(color_trait)
# rownames(color_trait) = species_names
# 
# # Trait tests
# tree_res = multi2di(tree)
# mps.result = multiPhylosignal(months, tree_res)
# 
# # Finding the monthly PD 
# monthly_phylo_dispersion = vector(mode="list", length=12)
# for(i in 1:12) {
#   month_species_list = rownames(months[months[i] == 1,])
#   month_community = community %>% select(which(colnames(community) %in% month_species_list))
#   month_tree = drop.tip(phy = tree, tip = setdiff(tree$tip.label,
#                                                       colnames(month_community)))
#   
#   ses.pd.result = ses.pd(month_community, month_tree, null.model="taxa.labels", runs=5)
#   monthly_phylo_dispersion[[i]] = ses.pd.result
#   print(paste0("Month ", i, " complete."))
# }
# 
# 
# # Finding the dominant color in each cell
# monthly_dominant_color_nums = c()
# for(i in 1:12) {
#   month_species_list = rownames(months[months[i] == 1,])
#   month_community = community %>% select(which(colnames(community) %in% month_species_list))
#   
#   color_table = data.table()
#   comm = t(data.frame(month_community))
#   rownames(comm) = colnames(month_community)
#   for(i in 1:nrow(month_community)) {
#     species = (comm[,i] == 1) %>% keep(comm[,i],.)
#     colors_cell = c(0,0,0,0,0,0,0,0,0,0)
#     if(length(species) > 1) {
#       for(j in 1:length(species)) {
#         species_name = names(species[j])
#         colors_species = unlist(colors[species_name,])
#         colors_cell = colors_cell + c(colors_species, 0)
#       }
#     } else { colors_cell = c(0,0,0,0,0,0,0,0,0,1) }
#     color_table = rbind(color_table, as.list(colors_cell), use.names=FALSE)
#     colnames(color_table) = c("red", "orange", "yellow", "green", "blue", "purple", 
#                               "pink", "white", "brown", "none")
#   }
#   dominant_colors = colnames(color_table)[max.col(color_table,ties.method="random")]
#   dom_cols_nums = c()
#   for(i in 1:length(dominant_colors)) {
#     dom_cols_nums[i] = which( colnames(color_table)==dominant_colors[[i]] )
#   }
#   monthly_dominant_color_nums = rbind(monthly_dominant_color_nums, dom_cols_nums)
# }
# 
# color_table = data.table()
# comm = t(data.frame(community))
# rownames(comm) = colnames(community)
# for(i in 1:nrow(community)) {
#   species = (comm[,i] == 1) %>% keep(comm[,i],.)
#   colors_cell = c(0,0,0,0,0,0,0,0,0,0)
#   if(length(species) > 0) {
#     for(j in 1:length(species)) {
#       species_name = names(species[j])
#       colors_species = unlist(colors[species_name,])
#       colors_cell = colors_cell + c(colors_species, 0)
#     }
#   } else { colors_cell = c(0,0,0,0,0,0,0,0,0,1) }
#   color_table = rbind(color_table, as.list(colors_cell), use.names=FALSE)
#   colnames(color_table) = c("red", "orange", "yellow", "green", "blue", "purple", 
#                             "pink", "white", "brown", "none")
# }
# dominant_colors = colnames(color_table)[max.col(color_table,ties.method="random")]
# dom_cols_nums = c()
# for(i in 1:length(dominant_colors)) {
#   dom_cols_nums[i] = which( colnames(color_table)==dominant_colors[[i]] )
# }
