# Create a Helical Plot for Time Series Data

Creates a helical plot to visualize time series data, emphasizing both
the magnitude and rate of change over time.

## Usage

``` r
heliPlot(
  data,
  group = NULL,
  facet = FALSE,
  xlabel = "Rate of Change",
  ylabel = "Values",
  arrow = TRUE,
  dateFont = 3,
  dateInterval = FALSE,
  sizeRange = c(1, 3),
  facetScales = "free",
  dateFormat = "%d %b %y",
  n = nrow(data),
  ...
)
```

## Arguments

- data:

  A data frame containing the time series data with required columns:
  "values_avg", "change_rate", and "date".

- group:

  (Optional) A string specifying the column name in \`data\` to use for
  grouping data in the plot. If NULL, no grouping is applied.

- facet:

  Logical indicating whether to facet the plot based on the \`group\`
  variable. If TRUE and \`group\` is NULL, an error is raised.

- xlabel:

  Label for the x-axis, defaults to "Rate of Change".

- ylabel:

  Label for the y-axis, defaults to "Values".

- arrow:

  Logical indicating whether to add an arrow to the end of each line,
  defaults to TRUE.

- dateFont:

  Numeric specifying the size of the date font, defaults to 3.

- dateInterval:

  Numeric specifying the interval at which date labels should be
  displayed. If FALSE, no date labels are shown.

- sizeRange:

  Numeric vector of length 2 specifying the range of line widths.

- facetScales:

  Character string indicating whether scales should be "fixed",
  "free_x", "free_y", or "free".

- dateFormat:

  Format for the date labels, defaults to d-b-y.

- n:

  Numeric specifying the number of points to interpolate along the
  spline, defaults to the number of rows in \`data\`.

- ...:

  Additional arguments passed on to \`ggplot2\` layer functions.

## Value

A \`ggplot\` object representing the helical plot.

## Examples

``` r
# Assuming `dataPrep` is a data frame prepared with the required structure:

if (FALSE) { # \dontrun{
heliPlot(dataPrep, group = "myGroup", arrow = TRUE, 
dateFont = 3, dateInterval = 30, sizeRange = c(1, 3))
} # }
```
