#' Prepare Data for Helical Plotting
#'
#' This function preprocesses a time series data for helical plotting by
#' applying a moving average and smoothing the values and their rate of change.
#' It also handles conversion of numeric dates to Date objects and ensures
#' proper alignment of the time series for plotting.
#'
#' @param dates A vector of dates associated with the values; can be numeric or
#'   Date objects. If numeric, they are treated as days since a given start date.
#' @param values A numeric vector of the time series values corresponding to the dates.
#' @param filterWidth The size of the moving window to calculate the moving average.
#'   Defaults to 7
#' @return A data frame suitable for helical plotting, containing the original
#'   dates, the smoothed values (`ch_avg`), the smoothed rate of change (`ch_rate`),
#'   and the endpoints for plotting (`yend`, `xend`).
#'
#' @examples
#' \dontrun{
#'   # Assume 'dates' and 'values' are available time series data
#'   prepared_data <- heliPrep(dates, values)
#'   # Now 'prepared_data' can be used for helical plotting with 'heliPlot'
#' }
#'
#' @export
#' @importFrom stats na.omit
heliPrep <- function(dates, values, filterWidth = 7) {

  # Validate inputs
  if (any(is.na(values))) {
    stop("time series contains internal NAs")
  }

  if (!is.numeric(values)) {
    stop("values must be numeric")
  }
  
  if (!(is.numeric(dates) || methods::is(dates,"Date"))) {
    stop("dates must be numeric or Date objects")
  }
  
  if (length(dates) != length(values)) {
    stop("dates and values must have the same length")
  }
  
  # Convert numeric dates to Date objects if necessary
  if (is.numeric(dates)) {
    start_date <- as.Date("1970-01-01")
    dates <- start_date + dates
  }

  m.factor <- as.integer((length(values)* filterWidth)/5)

  values.p <- c(rep(values[1], m.factor), values, rep(values[length(values)], m.factor))
  dates.p <- c(rep(dates[1], m.factor), dates, rep(dates[length(dates)], m.factor))

  # Calculate the moving average for values
  smooth_ma_values <- stats::filter(values.p, rep(1/filterWidth, filterWidth), sides = 2, method="convolution")
  
  # Calculate the rate of change for the moving average values
  ma_rate_change <- c(NA, diff(values.p))
  
  # Apply smoothing to the moving average values and rate of change
  smooth_ma_rate_change <- stats::filter(ma_rate_change, rep(1/filterWidth, filterWidth), sides = 2, method="convolution")

  # Derive the half moving window to cut dates
  halfWin <- ifelse(filterWidth%%2==1, (filterWidth-1)/2, filterWidth/2) 
  smooth_ma_values <- smooth_ma_values[m.factor:(length(smooth_ma_values) - m.factor)]
  smooth_ma_rate_change <- smooth_ma_rate_change[m.factor:(length(smooth_ma_rate_change) - m.factor)]
  valid_dates <- dates.p[m.factor:(length(dates.p) - (m.factor))]
      
  # Prepare the final data frame
  prepared_df <- data.frame(
    date = valid_dates,
    values_avg = smooth_ma_values,
    change_rate = smooth_ma_rate_change
    # yend = c(tail(smooth_ma_values, -1), NA),
    # xend = c(NA, head(smooth_ma_rate_change, -1))
  )
  
  # Remove rows with NAs
  prepared_df <- na.omit(prepared_df)
  return(prepared_df)
}