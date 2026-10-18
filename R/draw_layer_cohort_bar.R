#' Draw one or more numbered cohort bars
#' @param BT Data frame with `cohort`, `ymid`, `xlo`, and `xhi` columns.
#' @param pad Fraction of each cell dimension reserved on each edge for padding.
#' @param colour Outline, separator and text colour.
#' @param fill Cell background colour, or NA for transparency.
#' @param linewidth Width of outlines and separators.
#' @param group Optional column name identifying separate bars.
#' @param fontface Font face for cohort numbers, e.g. "plain" or "bold".
#' @param cellFill Optional column name in BT containing cell background colours.
#' NA entries retain the bar's background. Cell fills are below separators/text.
#' @param textColour Optional column name in BT containing per-cell text colours.
#' NA entries use the bar's text colour.
#' @return A list of ggplot2 layers, added to a plot with `+`.
#' @md
#' @export
draw_layer_cohort_bar <- function(BT, pad = 0.05, colour = "black",
                                  fill = NA, linewidth = 0.18, group = NULL,
                                  fontface = "plain", cellFill = NULL, textColour = NULL) {
  if (length(pad) != 1L || !is.finite(pad) || pad < 0 || pad >= 0.5) {
    stop("`pad` must be a fraction between 0 (inclusive) and 0.5 (exclusive).", call. = FALSE)
  }
  if (!is.numeric(linewidth) || length(linewidth) != 1L ||
      !is.finite(linewidth) || linewidth <= 0) {
    stop("`linewidth` must be one positive finite number.", call. = FALSE)
  }
  if (nrow(BT) == 0) return(list())
  for (column in list(cellFill, textColour)) {
    if (!is.null(column) && (length(column) != 1L || !is.character(column) ||
        is.na(column) || !column %in% names(BT))) {
      stop("Cell styling columns must name existing columns in BT.", call. = FALSE)
    }
  }
  groups <- if (is.null(group)) list(BT) else split(BT, BT[[group]])
  outlines <- do.call(rbind, lapply(groups, function(bar) {
    data.frame(xmin = min(bar$xlo), xmax = max(bar$xhi),
      ymin = min(bar$ymid) - 0.5, ymax = max(bar$ymid) + 0.5)
  }))
  backgrounds <- if (is.null(cellFill)) list() else {
    lapply(unique(stats::na.omit(BT[[cellFill]])), function(cell_colour) {
      cells <- BT[!is.na(BT[[cellFill]]) & BT[[cellFill]] == cell_colour, ]
      ggplot2::geom_rect(data = cells,
        ggplot2::aes(xmin = xlo, xmax = xhi, ymin = ymid - 0.5, ymax = ymid + 0.5),
        fill = cell_colour, colour = NA, inherit.aes = FALSE)
    })
  }
  label_colours <- if (is.null(textColour)) rep(colour, nrow(BT)) else BT[[textColour]]
  label_colours[is.na(label_colours)] <- colour
  labels <- lapply(unique(label_colours), function(label_colour) {
    ggfittext::geom_fit_text(data = BT[label_colours == label_colour, ],
      ggplot2::aes(xmin = xlo + pad * (xhi - xlo),
        xmax = xhi - pad * (xhi - xlo), ymin = ymid - 0.5 + pad,
        ymax = ymid + 0.5 - pad, label = cohort),
      colour = label_colour, fontface = fontface, min.size = 0, grow = TRUE, inherit.aes = FALSE,
      padding.x = grid::unit(0, "mm"), padding.y = grid::unit(0, "mm"))
  })
  list(
    ggplot2::geom_rect(data = outlines,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = fill, colour = colour, linewidth = linewidth, inherit.aes = FALSE),
    backgrounds,
    if (length(backgrounds)) ggplot2::geom_rect(data = outlines,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = colour, linewidth = linewidth, inherit.aes = FALSE) else NULL,
    ggplot2::geom_segment(data = BT,
      ggplot2::aes(x = xlo, xend = xhi, y = ymid - 0.5, yend = ymid - 0.5),
      colour = colour, linewidth = linewidth, inherit.aes = FALSE),
    labels
  )
}
