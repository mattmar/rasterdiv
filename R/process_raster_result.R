#' Process Raster Results
#'
#' This function processes the results of a list of calculations, packaging them into SpatRaster objects and naming them appropriately.
#'
#' @param out A list containing the results of calculations that need to be transformed into SpatRaster objects. Each element in the list corresponds to a different calculation result.
#' @param x A list of SpatRaster objects or a similar object used as a template for creating new SpatRaster objects. Specifically, `x[[1]]` is used as the template.
#' @param alpha Numeric or numeric vector indicating the weight(s) used in the distance calculations that generated 'out'. It is used to label the results appropriately.
#' @param window Numeric or numeric vector indicating the size(s) of the moving window(s) used in the calculations that generated 'out'. It is used to label the results.
#'
#' @return A list of SpatRaster objects corresponding to the different processed results. Each SpatRaster is named based on the 'alpha' and 'window' parameters used in the calculation. The naming convention is 'alpha.<alpha value>' for the inner lists and 'window.<window size>' for the outer list.
#'
#' @details The function is designed to post-process the results of spatial calculations performed on raster data. The typical use case is to process results from a function that performs calculations on different 'windows' of the data, using varying 'alpha' parameters, and returns the results as a list. This function takes that list, converts each element to a SpatRaster (using the first SpatRaster in 'x' as a template), and assigns appropriate names to each based on the 'alpha' and 'window' parameters.
#'
#' @examples
#' \dontrun{
#'
#' # Assume 'result_list' is obtained from a previous calculation, containing
#' # multiple results to be converted to SpatRaster objects.
#' # 'raster_template' is a list of SpatRaster objects used as templates.
#'
#' processed_results <- process_raster_result(out = result_list, 
#'                                            x = raster_template, 
#'                                            alpha = c(1, 2), 
#'                                            window = c(3, 5))
#' }
#'

process_raster_result <- function(out, x, alpha, window) {
  outR <- lapply(out, function(insm) {
    y <- lapply(insm, function(single_result) {
      terra::rast(single_result, type="", terra::crs(x[[1]]), terra::res(x[[1]]))
    })
    names(y) <- paste("alpha.", alpha, sep = "")
    return(y)
  })

  names(outR) <- paste("window.", window, sep = "")
  return(outR)
}
