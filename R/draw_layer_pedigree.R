#' Draw a positioned, numbered pedigree
#'
#' Draws direct child-to-parent lines, male squares and female circles.
#' Use the aspect ratio returned by the layout to preserve symbol shapes.
#' @param layout A list returned by [layout_pedigree_and_cohorts()].
#' @param colour Edge, node outline and label colour.
#' @param fill Node background colour.
#' @param linewidth Width of edges and outlines.
#' @param showGuides Draw thin dashed vertical guides through node centers,
#' behind the edges, filled nodes and labels.
#' @return A list of ggplot2 layers, added to a plot with `+`.
#' @md
#' @export
draw_layer_pedigree <- function(layout, colour = "black", fill = "white",
                                linewidth = 0.3, showGuides = TRUE) {
  nodes <- layout$nodes
  male <- nodes[nodes$sex == "M", ]
  female <- nodes[nodes$sex == "F", ]
  angle <- seq(0, 2 * pi, length.out = 65)
  circles <- if (nrow(female)) do.call(rbind, lapply(seq_len(nrow(female)), function(i) {
    data.frame(id = female$id[i],
      x = female$x[i] + (female$xmax[i] - female$xmin[i]) / 2 * cos(angle),
      y = female$y[i] + (female$ymax[i] - female$ymin[i]) / 2 * sin(angle))
  })) else data.frame(id = character(), x = numeric(), y = numeric())
  # Fit labels inside an inscribed box, including for circular nodes.
  nodes$text_xmin <- nodes$x + 0.68 * (nodes$xmin - nodes$x)
  nodes$text_xmax <- nodes$x + 0.68 * (nodes$xmax - nodes$x)
  nodes$text_ymin <- nodes$y + 0.68 * (nodes$ymin - nodes$y)
  nodes$text_ymax <- nodes$y + 0.68 * (nodes$ymax - nodes$y)
  list(
    if (showGuides) pedigree_guide_layer(layout, colour) else NULL,
    ggplot2::geom_segment(data = layout$edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = colour, linewidth = linewidth, inherit.aes = FALSE),
    ggplot2::geom_rect(data = male,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      colour = colour, fill = fill, linewidth = linewidth, inherit.aes = FALSE),
    ggplot2::geom_polygon(data = circles, ggplot2::aes(x = x, y = y, group = id),
      colour = colour, fill = fill, linewidth = linewidth, inherit.aes = FALSE),
    ggfittext::geom_fit_text(data = nodes,
      ggplot2::aes(xmin = text_xmin, xmax = text_xmax,
        ymin = text_ymin, ymax = text_ymax, label = label),
      colour = colour, min.size = 0, grow = TRUE, inherit.aes = FALSE,
      padding.x = grid::unit(0, "mm"), padding.y = grid::unit(0, "mm"))
  )
}

pedigree_guide_layer <- function(layout, colour) {
  ggplot2::geom_segment(data = layout$nodes,
    ggplot2::aes(x = x, xend = x), y = layout$ylim[1], yend = layout$ylim[2],
    colour = colour, linetype = "dashed", linewidth = 0.15, inherit.aes = FALSE)
}
