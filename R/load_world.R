#' Load Natural Earth world dataset
#'
#' This function loads and returns the World Vector data stored within the package.
#'
#' @return A `SpatVector` object representing the World vector data.
#' @details
#' SpatVector (EPSG: 4326) of the world dissolved on continents.
#' 
#' @examples 
#' \dontrun{
#' world_data <- load_world_vector()
#' plot(world_data)
#' }
#'
#' @source \code{https://www.naturalearthdata.com/}
#' 
#' @seealso \code{\link[terra]{vect}}
#' 
#' @export
load_world <- function() {
  # Set the path to the RDS file
  rds_file <- system.file("extdata", "world.rds", package = "rasterdiv")
  
  # Check if the file exists
  if (file.exists(rds_file)) {
    # Load the SpatVector object from the RDS file
    world_vector <- terra::vect(readRDS(rds_file))
  } else {
    stop("Data file 'world.rds' not found in the package.")
  }

  return(world_vector)
}
