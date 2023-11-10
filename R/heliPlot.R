#' Helical Plot for Time Series Data
#'
#' Creates a helical plot to visualize the relationship between rate of change and
#' cumulative averages in time series data. Can handle multiple time series by
#' faceting the plot based on a grouping variable.
#'
#' @param data A data frame containing the time series data with required columns: 'ch_avg', 'ch_rate', 'date', 'yend', and 'xend'.
#' @param group An optional string specifying the column name to use for grouping data in the plot.
#' @param facet Logical, if TRUE and a grouping variable is provided, will facet the plot based on the group.
#' @param xlabel Label for the x-axis.
#' @param ylabel Label for the y-axis.
#' @param arrow Logical, if TRUE, adds arrows to the plot lines to indicate direction.
#' @param textFont Numeric, the size of the text font for labels.
#' @param labelInterval Numeric, the interval at which labels are displayed. If FALSE, no labels are displayed.
#' @param sizeRange Numeric vector of size 2, indicating the range of line widths.
#' @param facetScales Character, if 'free', the scales of faceted plots are allowed to vary.
#' @param n Integer, the number of points to interpolate along the curve path.
#' @param ... Character, ggplot2 aesthetics (`aes()`).
#'
#' @return A `ggplot` object representing the helical plot with options for faceting and dynamic segment sizing.
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(
#'     date = seq(as.Date("2020-01-01"), by = "month", length.out = 12),
#'     ch_avg = runif(12, 5, 15),
#'     ch_rate = rnorm(12, 0, 1),
#'     yend = runif(12, 4, 16),
#'     xend = rnorm(12, -1, 1)
#'   )
#'   heliPlot(data, group = "date", facet = TRUE, xlabel = "Rate of Change",
#'            ylabel = "Values", arrow = TRUE, textFont = 3, labelInterval = 2,
#'            sizeRange = c(1, 3), facetScales = "free")
#' }
#' 
#' @export
#' @importFrom ggplot2 ggplot aes geom_segment geom_label scale_size scale_color_manual aes_string
#' @importFrom ggforce geom_link2
#' @importFrom viridis scale_colour_viridis
#' @importFrom ggplot2 theme_minimal theme element_text element_line unit ylab xlab
#' @importFrom ggplot2 geom_text geom_vline scale_colour_discrete scale_linewidth scale_colour_discrete
#' @importFrom ggplot2 facet_wrap guides 
#' @importFrom stats as.formula
#' @importFrom utils head tail
heliPlot <- function(data, group = NULL, facet = FALSE, xlabel = "Rate of Change", ylabel = "Values", arrow = TRUE, textFont = 3, labelInterval = FALSE, sizeRange = c(1, 3), facetScales="free", n=100, ...) {

  # Check for required columns
  required_cols <- c("ch_avg", "ch_rate", "date", "yend", "xend")
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
  p <- ggplot(data, aes(x = ch_rate, y = ch_avg), show.legend = FALSE) +
  group_aes +
  geom_vline(xintercept=0, col="grey", show.legend = FALSE) + 
  ggforce::geom_bspline2(
    aes(
      linewidth=ch_avg, 
      color=if(is.null(group)) {ch_rate} else get(group)),
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
      label_indices <- if(labelInterval) unlist(lapply(split(seq_along(data[[group]]), data[[group]]), function(x) {
        x[seq(1, length(x), by = labelInterval)]
        })) else 0
      arrow_indices <- unlist(lapply(split(seq_along(data[[group]]), data[[group]]), function(x) {
        x[c(1,length(x))]
        }))
      p <- if(labelInterval) p + 
      geom_text(data = data[label_indices, ],
       aes(label = format(date,"%d %b %y")), 
       size=textFont, vjust="inward",hjust="inward", colour="black", check_overlap=TRUE, show.legend = FALSE)
      # geom_point(data = data[label_indices, ], col="white",  show.legend = FALSE)
      else p
      } else {
        p <- if(labelInterval) p + 
        geom_text(data = data[seq(1, nrow(data), by = labelInterval), ],
         aes(label = format(date,"%b %y")), 
         size=textFont, vjust="inward",hjust="inward", colour="black", check_overlap=TRUE,
         show.legend = FALSE)
      # geom_point(data = data[seq(1, nrow(data), by = labelInterval), ], col="white", show.legend = FALSE)
      else p
      }

      if (facet & !is.null(group)) {
        p + 
        ggplot2::facet_wrap(as.formula(paste0("~", group)), scales = facetScales) + scale_colour_discrete()
        } else if (!facet & !is.null(group)) {
          p + scale_colour_discrete() 
         } else if (facet & is.null(group)) {
          stop("Faceting is true but no group variable is provided.")
          } else if ( is.null(group) ) {
            p + viridis::scale_colour_viridis(option="magma")
          }
        }