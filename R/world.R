#' Natural Earth world dataset
#'
#' A \code{SpatVector} (EPSG: 4326) of the continents.
#'
#' This dataset represents the world, with spatial information dissolved on continents.
#' It is stored as a \code{SpatVector}, suitable for various spatial operations and mapping.
#'
#' @format A \code{SpatVector} containing the following columns:
#' \describe{
#'   \item{world}{SpatVector of the world dissolved on continents. Details about columns should be listed here if applicable.}
#' }
#' @usage load_world()
#' @source \url{https://www.naturalearthdata.com/}
#' @references \url{https://www.naturalearthdata.com/}
#' @name world
#' @docType data
#' @keywords datasets
#' @examples
#' world <- readRDS(system.file("extdata", "world.rds", package = "rasterdiv"))
NULL