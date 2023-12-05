#' Create a Helical Plot for Time Series Data
#'
#' @description
#' Creates a helical plot to visualize time series data, emphasizing both
#' the magnitude and rate of change over time.
#'
#' @param data A data frame containing the time series data with required columns: "values_avg", "change_rate", and "date".
#' @param group (Optional) A string specifying the column name in `data` to use for grouping data in the plot. If NULL, no grouping is applied.
#' @param facet Logical indicating whether to facet the plot based on the `group` variable. If TRUE and `group` is NULL, an error is raised.
#' @param xlabel Label for the x-axis, defaults to "Rate of Change".
#' @param ylabel Label for the y-axis, defaults to "Values".
#' @param arrow Logical indicating whether to add an arrow to the end of each line, defaults to TRUE.
#' @param dateFont Numeric specifying the size of the date font, defaults to 3.
#' @param dateInterval Numeric specifying the interval at which date labels should be displayed. If FALSE, no date labels are shown.
#' @param sizeRange Numeric vector of length 2 specifying the range of line widths.
#' @param facetScales Character string indicating whether scales should be "fixed", "free_x", "free_y", or "free".
#' @param dateFormat Format for the date labels, defaults to d-b-y.
#' @param n Numeric specifying the number of points to interpolate along the spline, defaults to the number of rows in `data`.
#' @param ... Additional arguments passed on to `ggplot2` layer functions.
#'
#' @return A `ggplot` object representing the helical plot.
#'
#' @examples
#' # Assuming `dataPrep` is a data frame prepared with the required structure:
#'
#' \dontrun{
#' heliPlot(dataPrep, group = "myGroup", arrow = TRUE, 
#' dateFont = 3, dateInterval = 30, sizeRange = c(1, 3))
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom ggforce geom_bspline2
#' @importFrom viridis scale_colour_viridis

heliPlot <- function(data, group = NULL, facet = FALSE, xlabel = "Rate of Change", ylabel = "Values", arrow = TRUE, dateFont = 3, dateInterval = FALSE, sizeRange = c(1, 3), facetScales="free", dateFormat="%d %b %y", n=nrow(data), ...) {

  # Check for required columns
  required_cols <- c("values_avg", "change_rate", "date")
  if (!all(required_cols %in% colnames(data))) {
    stop("Data frame lacks required columns.")
  }

  # Check if group variable is provided and is in data
  if (!is.null(group) && !group %in% colnames(data)) {
    stop("Grouping variable not found in data.")
  }

  # Set group aesthetic if group is provided
  group_aes <- if (is.null(group)) aes_string(group=1, ...) else aes_string(group = group, ...) 

  # Create plot
  p <- ggplot(data, aes(x = data[,3], y = data[,2]), show.legend = FALSE) +
  group_aes +
  geom_vline(xintercept=0, col="grey", show.legend = FALSE) + 
  ggforce::geom_bspline2(
    aes(
    linewidth=data[,2], 
    color=if(is.null(group)) {data[,3]} else get(group)),
    alpha=0.7,
    show.legend = TRUE, 
    lineend = 'round', 
    n = n, 
    arrow = if(arrow) arrow(type = 'closed', ends = "last", length = unit(0.05, "inches")) else NULL) +
  scale_linewidth(range = sizeRange) +
  ylab(ylabel) + 
  xlab(xlabel) +
  theme_minimal(base_size = 14) +
  guides(linewidth="none") +
  theme(
    legend.position = "none",
    axis.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(size = 9),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.text.x = ggplot2::element_text(size = 15),
    panel.grid.major.x = ggplot2::element_line(color = "grey80", size = 0.5,linetype=2),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    aspect.ratio = 1)

    # If group is not NULL, calculate label indices based on group
    if (!is.null(group)) {
      label_indices <- if(dateInterval) unlist(lapply(split(seq_along(data[[group]]), data[[group]]), function(x) {
        x[seq(1, length(x), by = dateInterval)]
        })) else 0
      arrow_indices <- unlist(lapply(split(seq_along(data[[group]]), data[[group]]), function(x) {
        x[c(1,length(x))]
        }))
      p <- if(dateInterval) p + 
      geom_text(data = data[label_indices, ],
       aes(label = format(date,dateFormat)), 
       size=dateFont, vjust="inward",hjust="inward", colour="black", check_overlap=TRUE, show.legend = FALSE)
      else p
      } else {
        p <- if(dateInterval) p + 
        geom_text(data = data[seq(1, nrow(data), by = dateInterval), ],
         aes(label = format(date,dateFormat)), 
         size=dateFont, vjust="inward",hjust="inward", colour="black", check_overlap=TRUE,
         show.legend = FALSE)
      # geom_point(data = data[seq(1, nrow(data), by = dateInterval), ], col="white", show.legend = FALSE)
      else p
      }

      if (facet & !is.null(group)) {
        p + 
        ggplot2::facet_wrap(stats::as.formula(paste0("~", group)), scales = facetScales) + scale_colour_discrete()
        } else if (!facet & !is.null(group)) {
          p + scale_colour_discrete() 
         } else if (facet & is.null(group)) {
          stop("Faceting is true but no group variable is provided.")
          } else if ( is.null(group) ) {
            p + viridis::scale_colour_viridis(option="magma")
          }
        }