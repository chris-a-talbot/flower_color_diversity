# TO DO: Standardize so all output has only the value of interest for null_distributions
get_cell_distributions = function(null_distributions) {
  for(i in 1:length(null_distributions)) {
    # Temporary list to store the FunRao values for the current overarching item
    tempResults <- list()
    
    cell_num = length(null_distributions[[i]][[1]][[1]])
    
    # Initialize lists to collect FunRao values for each label
    for(j in 1:cell_num) {
      tempResults[[as.character(j)]] <- vector()
    }
    
    # Iterate over each of the 100 items within the current overarching item
    for(item in null_distributions[[i]]) {
      # Access the FunRao list
      val_list <- item$Q
      
      # Iterate over each integer label in the FunRao list
      for(label in 1:length(val_list)) {
        # Append the value to the appropriate list in tempResults
        tempResults[[as.character(label)]] <- c(tempResults[[as.character(label)]], val_list[[label]])
      }
    }
    
    # Add the compiled FunRao values for the current overarching item to the results list
    cell_distributions_1[[as.character(i)]] <- tempResults
  }
}

# REQUIRES: A cell_distributions object, with each of n cells being a list of values or NA
# REQUIRES: A number of cores >= 1
# MODIFIES: N/A
# EFFECTS: Outputs a list of ecdf functions, one for each distribution in cell_distributions
get_ecdf_list = function(cell_distributions, cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  if(cores == 1) {
    ecdf_list = lapply(1:length(cell_distributions), function(cell_num) {
      row_data <- unlist(cell_distributions[cell_num])
      if(all(is.na(row_data))) {
        return(NA)
      } else {
        # Remove NAs and compute ECDF, then add to the list
        return(ecdf(na.omit(row_data)))
      }
    })
  } else {
    ecdf_list = mclapply(1:length(cell_distributions), function(cell_num) {
      row_data <- unlist(cell_distributions[cell_num])
      if(all(is.na(row_data))) {
        return(NA)
      } else {
        # Remove NAs and compute ECDF, then add to the list
        return(ecdf(na.omit(row_data)))
      }
    }, mc.cores=cores)
  }
  
  # Retain site index info
  names(ecdf_list) = names(cell_distributions)
  
  return(ecdf_list)
}

# REQUIRES: A cell_distributions_list object, which is a list of lists of 
#           distributions for each cell
# REQUIRES: A number of cores >= 1
# MODIFIES: N/A
# EFFECTS: Outputs a list of lists of ecdf functions, one for each distibution 
#          in each list of cell_distributions
get_ecdf_lists = function(cell_distributions_list, cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  if(cores == 1) {
    ecdfs_list = lapply(1:length(cell_distributions_list), 
                        function(cell_num) get_ecdf_list(cell_distributions_list[[cell_num]]))
  } else {
    ecdfs_list = mclapply(1:length(cell_distributions_list),
                          function(cell_num) get_ecdf_list(cell_distributions_list[[cell_num]], 
                                                           cores=cores), 
                          mc.cores=cores)
  }
  
  return(ecdfs_list)
}

# REQUIRES: A list of ecdf functions
# REQUIRES: A list of observed values, of equal length
# MODIFIES: N/A
# EFFECTS: For each ecdf function, plugs in corresponding observed value - ECDF(observed_val)
# EFFECTS: Outputs a list of ECDF values for each cell
get_ecdf_vals = function(ecdf_list, observed_vals) {
  # Gets a list of ecdfs for each cell for a given month
  output_list <- vector("list", length = length(ecdf_list))
  names(output_list) <- names(ecdf_list)
  
  output_list = lapply(1:length(ecdf_list), function(ecdf_index) {
    if (is.na(ecdf_list[[ecdf_index]])) {
      # If the ecdf function is NA, set the output to NA
      return(NA)
    } else {
      cell_name <- names(ecdf_list)[ecdf_index]
      if (!cell_name %in% names(observed_vals)) {
        # If there's no observed value for this cell, set the output to NA
        return(NA)
      } else {
        # Calculate the difference between observed value and the ecdf function range
        observed <- observed_vals[[cell_name]]
        ecdf_func <- ecdf_list[[ecdf_index]]
        
        # Get the probabiliy of the observed value being at least as large as it is
        p = ecdf_func(observed)
        
        # Add to list
        return(p)
      }
    }
  })
  
  names(output_list) = names(ecdf_list)
  return(output_list)
}

# REQUIRES: A resolution for the raster
# REQUIRES: A location of the .tif file
# REQUIRES: A month to determine which .tif to select
# MODIFIES: N/A
# EFFECTS: Outputs a list c(nrow, ncol) of the dimensions of the raster
get_dims_from_tif = function(res=1, location="./data", month="jul") {
  raster = rast(paste0(location, "/", month, "_stack_", as.character(res), ".tif"))
  return(c(nrow(raster), ncol(raster)))
}

# REQUIRES: A list of ecdf_vals (or really any values of equal length to the raster being made)
# REQUIRES: Dims = c(nrow, ncol) for the raster; if NA, tries to get dims using month
# REQUIRES: A resolution for the raster
# REQUIRES: An extent for the raster
# REQUIRES: A month that the raster is representing, used to get dims; if NA, uses July as a default
# REQUIRES: A plot label
# MODIFIES: N/A
# EFFECTS: Outputs a plot of the input values on a raster of the given size, laid over a world map
plot_ecdf = function(ecdf_vals, dims=NA, res=1, ext=c(-92.5, -67.5, 38, 48), month="jul", 
                     plot_label="Flowers") {
  
  map = ne_countries(scale='medium', type='map_units', returnclass='sf', 
                     continent=c("North America")) # Generic world map
  
  if(is.na(dims)) {
    if(is.na(month)) {
      stop("No dimensions supplied for plot_ecdf()!")
    } else {
      dims = get_dims_from_tif(res=res, month=month)
    }
  }

  raster = rast(nrow=dims[1], ncol=dims[2], xmin=ext[1], xmax=ext[2], 
                ymin=ext[3], ymax=ext[4])
  values(raster) = unlist(ecdf_vals)
  
  # Convert SpatRaster to data frame for ggplot2
  raster_df <- as.data.frame(terra::as.data.frame(raster, xy=TRUE), na.rm=TRUE)
  names(raster_df) <- c("x", "y", "value")
  p <- ggplot() + geom_sf(data=map, color="black") + 
    coord_sf(xlim=c(-92.5, -67.5), ylim=c(38, 48)) +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value), alpha=0.75) +
    scale_fill_gradient2(low = "red", mid="yellow", high = "green", name = "ECDF(FRic)", 
                         limits = c(0, 1), breaks = c(0, 1), labels = c("0", "1"), midpoint=0.5) +
    theme_minimal() +
    labs(title = plot_label, subtitle="") +
    theme(plot.title = element_text(hjust = 0.5),# Center the title
          plot.subtitle = element_text(hjust=0.5),
          axis.text.x = element_blank(), # Remove x-axis text
          axis.ticks.x = element_blank(), # Optionally, remove x-axis ticks if desired
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title.x.bottom = element_blank(),
          axis.title.y.left = element_blank()) # Adjust x-axis title position if necessary
  
  return(p)
}


