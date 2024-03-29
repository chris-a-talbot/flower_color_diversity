# REQUIRES: A presence/absence matrix of species across sites
# REQUIRES: A list of site indices to calculate from
# MODIFIES: N/A
# EFFECTS: Returns the species richness of each site
get_species_richness = function(pa_matrix, indices=NULL) {
  
  # Use all sites, if no sites given
  if(!is.null(indices)) {
    pa_matrix = pa_matrix[indices,]
  }
  
  # Get the species richness of each site
  richness = rowSums(pa_matrix)
  
  # Set the names of each item to reflect the site index
  if(!is.null(indices)) {
    names(richness) = indices
  } else {
    names(richness) = 1:length(richness)
  }
  
  return(richness)
}

# REQUIRES: A presence/absence matrix of species across sites
# REQUIRES: A list of site indices to calculate from
# REQUIRES: A data.table of color categories for each species
# REQUIRES: A logical value determining whether to use a random assemblage of species
# REQUIRES: A number of cores to use for parallel processing
# REQUIRES: A seed for randomly selecting species for random assemblages
# MODIFIES: Seed
# EFFECTS: Returns the sum of color categories in each site
get_color_distributions = function(pa_matrix, indices=NULL, colors, null_hypothesis=FALSE, cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  # Use all sites if no indices provided
  if(is.null(indices)) {
    indices = 1:nrow(pa_matrix)
  }
  
  # Get the summed trait lists for every site
  summed_traits_list = mclapply(indices, function(i) {
    # Identify the species present at this site
    species_present = colnames(pa_matrix[,which(pa_matrix[i, ] == 1)])
    
    # Filter the trait data to the species of interest, or a random assemblage
    if(null_hypothesis) {
      filtered_traits = colors[sample(nrow(colors), length(species_present)),]
    } else {
      filtered_traits = colors[species_present,]
    }
    
    # Sum up the trait values column-wise
    summed_traits = colSums(filtered_traits)
    
    return(summed_traits)
  }, mc.cores=cores)
  
  summed_trait_df = do.call(rbind, mclapply(summed_traits_list, function(df) {
    return(as.data.frame(t(unlist(df))))
  }, mc.cores=cores))
  
  rownames(summed_trait_df) = indices
  
  return(summed_trait_df)
}

# REQUIRES: A list of color distributions from get_color_distributions()
# REQUIRES: A number of cores to use for parallel processing
# MODIFIES: N/A
# EFFECTS: Returns the Simpson's Evenness for each site
get_simpsons_evenness = function(color_distributions, cores=1) {
  
  # Calculate Simpson's Evenness for each site
  simpsons_e = mclapply(1:nrow(color_distributions), function(i) {
    return(simpson_e(color_distributions[i,]))
  }, mc.cores=cores)
  
  # Reformat the data
  simpsons_e = unlist(simpsons_e)
  
  # Retain site indices
  names(simpsons_e) = rownames(color_distributions)
  
  return(simpsons_e)
}
