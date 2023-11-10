#' Load Copernicus NDVI data
#'
#' This function loads and returns the Copernicus Long Term (1999-2017) NDVI Overview stored within the package.
#'
#' @return A `SpatRaster` object representing the Copernicus NDVI data.
#' @export
load_copNDVI <- function() {
  # Set the path to the RDS file
  rds_file <- system.file("extdata", "copNDVI.rds", package = "rasterdiv")
  
  # Check if the file exists
  if (file.exists(rds_file)) {
    # Load the SpatRaster object from the RDS file
    copNDVI <- terra::rast(readRDS(rds_file))
  } else {
    stop("Data file 'copNDVI.rds' not found in the package.")
  }
  
  return(copNDVI)
}
