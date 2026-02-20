# Prepare Data for Helical Plotting

This function preprocesses a time series data for helical plotting by
applying a moving average and smoothing the values and their rate of
change. It also handles conversion of numeric dates to Date objects and
ensures proper alignment of the time series for plotting.

## Usage

``` r
heliPrep(dates, values, filterWidth = 7)
```

## Arguments

- dates:

  A vector of dates associated with the values; can be numeric or Date
  objects. If numeric, they are treated as days since a given start
  date.

- values:

  A numeric vector of the time series values corresponding to the dates.

- filterWidth:

  The size of the moving window to calculate the moving average.
  Defaults to 7

## Value

A data frame suitable for helical plotting, containing the original
dates, the smoothed values (\`ch_avg\`), the smoothed rate of change
(\`ch_rate\`), and the endpoints for plotting (\`yend\`, \`xend\`).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Assume 'dates' and 'values' are available time series data
  prepared_data <- heliPrep(dates, values)
  # Now 'prepared_data' can be used for helical plotting with 'heliPlot'
} # }
```
