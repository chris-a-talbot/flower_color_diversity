################################################################################
##       This is a source file for raster-based helper functions.             ##
##                      Chris Talbot, March 2024                              ##
################################################################################

# REQUIRES: A directory (relative to working) that contains shapefiles
# MODIFIES: N/A
# EFFECTS: Returns a list of all valid shapefile file names in the given directory
get_shapefile_names = function(location="./shp") {
  shapefile_names = list.files(path=location, pattern="\\.shp$")
  return(shapefile_names)
}

# REQUIRES: A list of shapefile file names
# REQUIRES: The directory (relative to working) containing the shapefiles
# REQUIRES: A number of cores >=1 - runs on mclapply if cores>1, requires Unix
# REQUIRES: A debug logical value
# MODIFIES: N/A
# EFFECTS: Returns a list of all shapefiles
# EFFECTS: If debug==TRUE, prints to console each time a shapefile is read
get_species_ranges = function(file_names, location="./shp", cores=1, debug=FALSE) {
  # cores must be >= 1
  if(cores < 1) { 
    stop("Invalid core number for get_species_ranges - must be >1")
  }
  
  # REQUIRES: A file path of a valid shapefile
  # MODIFIES: N/A
  # EFFECTS: Flattens all regions in the sf object read in into a single region
  # EFFECTS: Returns the flattened region as an sf object 
  get_flat_sf_from_file = function(file_path) {
    # Read the shapefile
    sf_obj <- read_sf(paste0(location, "/", file_path))
    
    # Dissolve all geometries into a single geometry for this shapefile
    dissolved_geom <- st_union(sf_obj$geometry)
    
    # Create a new sf object with the dissolved geometry
    dissolved_sf <- st_sf(geometry = st_sfc(dissolved_geom))
    
    if(debug) {
      print(paste(file_path, "shapefile read."))
    }
    
    return(dissolved_sf)
  }
  
  # Run on lapply if running on 1 core
  if(cores == 1) {
    species_ranges = lapply(file_names, 
                            function(file_name) get_flat_sf_from_file(file_name))
    names(species_ranges) = file_path_sans_ext(file_names)
  # Run on mclapply if running on >1 core
  } else {
    species_ranges = mclapply(file_names, 
                              function(file_name) get_flat_sf_from_file(file_name),
                              mc.cores=cores)
    names(species_ranges) = file_path_sans_ext(file_names)
  }
  
  return(species_ranges)
}

# REQUIRES: A list of species occurrences in a single data.table with columns "species", "decimalLatitude", and "decimalLongitude"
# REQUIRES: A data.frame of phenological traits for all species of interest,
#           with a "species" column and a column for every month
# REQUIRES: A desired extent for the rasters, in the form c(ymin, ymax, xmin, xmax)
# REQUIRES: A desired resolution for the rasters, default=1
# REQUIRES: A logical value for write, default=FALSE
# REQUIRES: A location to write files to if write=TRUE
# REQUIRES: A list of months to run through, default=feb:nov
# REQUIRES: A number of cores >=1 - runs on mclapply if cores>1, requires Unix
# MODIFIES: N/A
# EFFECTS: Returns a list of rasters for each month, with each raster having a 
#          layer for each species in that month on the phenological traits data.frame.
#          These rasters do NOT work with parallel computing or saving RData -
#          use get_monthly_rasters_from_tif() and convert to p/a matrix for this.
# EFFECTS: If write=TRUE, writes each monthly raster to the specified location
get_monthly_rasters_from_occurrences = function(occs, pheno, ext, res=0.09, write=FALSE, 
                                           overwrite=TRUE, location="./data",
                                           months=c("feb", "mar", "apr", "may", 
                                                    "jun", "jul", "aug", "sep", 
                                                    "oct", "nov"), 
                                           cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  # REQUIRES: A 3-letter abbreviated month
  # MODIFIES: N/A
  # EFFECTS: Returns a raster for a given month
  get_raster_for_month = function(month) {
    month_species_list = pheno %>% filter(!!as.symbol(month) == 1) %>% 
      pull(species)
    month_species_occs = occs[species %in% month_species_list]
    
    if(length(month_species_occs) > 0) { 
      
      # Convert 'decimalLatitude' and 'decimalLongitude' to 'Y' and 'X' respectively, for each species
      species_list_now = month_species_occs[, .(X = decimalLongitude, Y = decimalLatitude), by = species]
      
      # Convert to a list of data.tables, one per species
      species_data_tables_list = split(species_list_now, species_list_now$species)
      
      # Remove the 'species' column from each data.table in the list
      species_data_tables_list = lapply(species_data_tables_list, function(dt) dt[, species := NULL])
      
      names(species_data_tables_list) = month_species_list
      
      stack <- rast(rasterStackFromOccurrences(species_data_tables_list, resolution=res, extent=ext))
      
      names(stack) = month_species_list
    } else { 
      stop("Invalid data given to get_raster_for_month() - 0 species ranges found")
    }
    
    if(write) {
      terra::writeRaster(stack, paste0(location, "/", month, "_stack_", as.character(res), ".tif"), 
                  overwrite=overwrite)
    }
    
    return(stack)
  }
  
  # Run on lapply if running on 1 core
  if(cores == 1) {
    monthly_rasters = lapply(months, function(month) get_raster_for_month(month))
    # Run on mclapply if running on >1 core
  } else {
    monthly_rasters = mclapply(months, function(month) get_raster_for_month(month),
                               mc.cores=cores)
  }
  
  return(monthly_rasters)
}

# REQUIRES: A list of species ranges as read-in shapefiles (output of get_species_ranges)
# REQUIRES: A data.frame of phenological traits for all species of interest,
#           with a "species" column and a column for every month
# REQUIRES: A desired extent for the rasters, in the form c(ymin, ymax, xmin, xmax)
# REQUIRES: A desired resolution for the rasters, default=1
# REQUIRES: A logical value for write, default=FALSE
# REQUIRES: A location to write files to if write=TRUE
# REQUIRES: A list of months to run through, default=feb:nov
# REQUIRES: A number of cores >=1 - runs on mclapply if cores>1, requires Unix
# MODIFIES: N/A
# EFFECTS: Returns a list of rasters for each month, with each raster having a 
#          layer for each species in that month on the phenological traits data.frame.
#          These rasters do NOT work with parallel computing or saving RData -
#          use get_monthly_rasters_from_tif() and convert to p/a matrix for this.
# EFFECTS: If write=TRUE, writes each monthly raster to the specified location
get_monthly_rasters_from_ranges = function(ranges, pheno, ext, res=1, write=FALSE, 
                                           overwrite=TRUE, location="./data",
                                           months=c("feb", "mar", "apr", "may", 
                                                    "jun", "jul", "aug", "sep", 
                                                    "oct", "nov"), 
                                           cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  # REQUIRES: A 3-letter abbreviated month
  # MODIFIES: N/A
  # EFFECTS: Returns a raster for a given month
  get_raster_for_month = function(month) {
    month_species_list = pheno %>% filter(!!as.symbol(month) == 1) %>% 
      pull(species)
    month_species_ranges = ranges[names(ranges) %in% month_species_list]
    
    if(length(month_species_ranges) > 0) { 
      stack <- rasterStackFromPolyList(month_species_ranges, 
                                       resolution=res, 
                                       retainSmallRanges=FALSE, 
                                       extent=ext) 
      if(length(month_species_ranges) == 1) { 
        names(stack) <- month_species_list[1] 
      }
    } else { 
      stop("Invalid data given to get_raster_for_month() - 0 species ranges found")
    }
    
    if(write) {
      writeRaster(stack, paste0(location, "/", month, "_stack_", as.character(res), ".tif"), 
                  overwrite=overwrite)
    }
    
    return(stack)
  }
  
  # Run on lapply if running on 1 core
  if(cores == 1) {
    monthly_rasters = lapply(months, function(month) get_raster_for_month(month))
  # Run on mclapply if running on >1 core
  } else {
    monthly_rasters = mclapply(months, function(month) get_raster_for_month(month),
                               mc.cores=cores)
  }
  
  return(monthly_rasters)
}

# REQUIRES: A raster containing >= 1 layer representing species across sites
# MODIFIES: N/A
# EFFECTS: Outputs a presence/absence matrix for each layer across all sites
raster_to_pa_matrix <- function(raster_stack) {
  # Use the values function from terra package to extract raster values
  raster_values <- values(raster_stack)
  
  # Convert to binary presence-absence data
  pa_matrix <- as.data.frame(ifelse(!is.na(raster_values) & raster_values > 0, 1, 0))
  
  # Set column names to species names (assuming layer names in the stack are species names)
  colnames(pa_matrix) <- names(raster_stack)
  
  return(pa_matrix)
}

# REQUIRES: A 3-letter abbreviated month
# REQUIRES: A resolution of the raster
# REQUIRES: A location containing the raster
# MODIFIES: N/A
# EFFECTS: Outputs a SpatRaster object for the given month/resolution combo
get_month_raster_from_tif = function(months, res=1, location="./data", cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  if(cores == 1) {
    monthly_rasters = lapply(months, 
                             function(month) return(rast(paste0(location, "/", 
                                                                month, "_stack_", 
                                                                as.character(res), 
                                                                ".tif"))))
  } else {
    monthly_rasters = mclapply(months, 
                             function(month) return(rast(paste0(location, "/", 
                                                                month, "_stack_", 
                                                                as.character(res), 
                                                                ".tif"))),
                             mc.cores=cores)
  }
  
  names(monthly_rasters) = months
  
  if(length(monthly_rasters) == 1) { 
    monthly_rasters = monthly_rasters[[1]]
  } 
  
  return(monthly_rasters)
}

# REQUIRES: A 3-letter abbreviated month
# REQUIRES: A resolution of the raster
# REQUIRES: A location containing the raster
# MODIFIES: N/A
# EFFECTS: Outputs a data.frame presence/absence matrix for the given month/resolution combo
get_month_pa_matrix_from_tif = function(months, res=1, location="./data", cores=1) {
  # Cores must be >= 1
  if(cores < 1) {
    stop("Invalid core number for get_monthly_rasters - must be >1")
  } 
  
  if(cores==1) {
    monthly_pa_matrices = lapply(months, function(month) {
      matrix = raster_to_pa_matrix(get_month_raster_from_tif(month, res=res, 
                                                             location=location))
      return(matrix)
    })
  } else {
    monthly_pa_matrices = mclapply(months, function(month) {
      matrix = raster_to_pa_matrix(get_month_raster_from_tif(month, res=res, 
                                                             location=location))
      return(matrix)
    }, mc.cores=cores) 
  }
  
  names(monthly_pa_matrices) = months
  
  if(length(monthly_pa_matrices) == 1) { 
    monthly_pa_matrices = monthly_pa_matrices[[1]]
  } 
  
  return(monthly_pa_matrices)
}
