# Plotting CKMR model graphs ---------------------------------------------------
#
# The plotting pipeline has two explicit steps:
#   1. layout_ckmr_model() turns a graph-like CKMR model into coordinates.
#   2. plot_ckmr_model() renders those coordinates with ordinary ggplot layers.
#
# Keeping these separate should make it easier to add alternate layouts, Shiny
# interactions, path highlighting, and eventually animation without changing the
# core model object.

layout_ckmr_model <- function(model, layout = "two_stage_po") {
  stopifnot_ckmr_model(model)

  if (layout != "two_stage_po") {
    stop("Only `layout = \"two_stage_po\"` is implemented so far.", call. = FALSE)
  }

  nodes <- model$nodes
  edges <- model$edges

  # The first layout is intentionally hand-authored and pedagogical. Later,
  # layout functions can become more algorithmic for time-by-age grids.
  coords <- data.frame(
    node_id = c(
      "adult_F", "adult_M", "mating", "offspring_pool",
      "adult_F_sample", "adult_M_sample", "offspring_sample"
    ),
    x = c(0, 4, 2, 2, -1.10, 2.90, 0.78),
    y = c(3, 3, 1.85, 0.35, 3.68, 3.68, 1.06),
    width = c(1.25, 1.25, 0.80, 1.35, 0.28, 0.28, 0.28),
    height = c(1.25, 1.25, 0.36, 0.95, 0.28, 0.28, 0.28),
    stringsAsFactors = FALSE
  )

  nodes <- merge(nodes, coords, by = "node_id", all.x = TRUE, sort = FALSE)
  nodes <- nodes[match(model$nodes$node_id, nodes$node_id), ]

  missing_coords <- nodes$node_id[is.na(nodes$x)]
  if (length(missing_coords) > 0) {
    stop(
      "The two-stage PO layout has no coordinates for: ",
      paste(missing_coords, collapse = ", "),
      call. = FALSE
    )
  }

  edges <- merge(
    edges,
    nodes[, c("node_id", "x", "y", "width", "height")],
    by.x = "from",
    by.y = "node_id",
    all.x = TRUE,
    sort = FALSE
  )
  names(edges)[names(edges) %in% c("x", "y", "width", "height")] <-
    c("x_from", "y_from", "width_from", "height_from")

  edges <- merge(
    edges,
    nodes[, c("node_id", "x", "y", "width", "height")],
    by.x = "to",
    by.y = "node_id",
    all.x = TRUE,
    sort = FALSE
  )
  names(edges)[names(edges) %in% c("x", "y", "width", "height")] <-
    c("x_to", "y_to", "width_to", "height_to")

  # Anchor edges on the visible boundary of the node shapes. This is already
  # enough to make the grammar feel diagrammatic rather than like a generic
  # node-link graph.
  edges <- add_edge_anchors(edges)

  list(
    model = model,
    nodes = nodes,
    edges = edges,
    offspring_shapes = make_offspring_shapes(nodes),
    mini_panels = make_mini_panels(model, nodes)
  )
}

plot_ckmr_model <- function(model,
                            overlay = c("none", "abundance", "fecundity", "length_distribution"),
                            layout = "two_stage_po") {
  overlay <- match.arg(overlay)
  arranged <- layout_ckmr_model(model, layout = layout)

  nodes <- arranged$nodes
  edges <- arranged$edges

  population_nodes <- nodes[nodes$node_type == "population", ]
  reproduction_nodes <- nodes[nodes$node_type == "reproduction", ]
  sample_nodes <- nodes[nodes$node_type == "sample", ]

  ggplot2::ggplot() +
    draw_reproduction_edges(edges) +
    draw_sampling_edges(edges) +
    draw_population_nodes(population_nodes, overlay) +
    draw_reproduction_nodes(reproduction_nodes) +
    draw_offspring_pool(arranged$offspring_shapes) +
    draw_sample_nodes(sample_nodes) +
    draw_mini_panels(arranged$mini_panels, overlay) +
    ggplot2::coord_equal(xlim = c(-1.25, 5.0), ylim = c(-0.45, 4.25), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "none"
    )
}

add_edge_anchors <- function(edges) {
  edges$x_start <- edges$x_from
  edges$y_start <- edges$y_from
  edges$x_end <- edges$x_to
  edges$y_end <- edges$y_to

  repro_input <- edges$edge_type == "reproduction_input"
  birth_output <- edges$edge_type == "birth_output"
  sampling <- edges$edge_type == "sampling"

  # Reproduction inputs leave from the lower inside corners of the adult boxes
  # and enter the upper corners of the mating black box.
  left_input <- repro_input & edges$x_from < edges$x_to
  right_input <- repro_input & edges$x_from > edges$x_to

  edges$x_start[left_input] <- edges$x_from[left_input] + edges$width_from[left_input] / 2
  edges$y_start[left_input] <- edges$y_from[left_input] - edges$height_from[left_input] / 2
  edges$x_end[left_input] <- edges$x_to[left_input] - edges$width_to[left_input] / 2
  edges$y_end[left_input] <- edges$y_to[left_input] + edges$height_to[left_input] / 2

  edges$x_start[right_input] <- edges$x_from[right_input] - edges$width_from[right_input] / 2
  edges$y_start[right_input] <- edges$y_from[right_input] - edges$height_from[right_input] / 2
  edges$x_end[right_input] <- edges$x_to[right_input] + edges$width_to[right_input] / 2
  edges$y_end[right_input] <- edges$y_to[right_input] + edges$height_to[right_input] / 2

  # Birth output drops from the mating box into the offspring pot.
  edges$x_start[birth_output] <- edges$x_from[birth_output]
  edges$y_start[birth_output] <- edges$y_from[birth_output] - edges$height_from[birth_output] / 2
  edges$x_end[birth_output] <- edges$x_to[birth_output]
  edges$y_end[birth_output] <- edges$y_to[birth_output] + edges$height_to[birth_output] / 2

  # Samples emerge from the northwest side of whatever population component is
  # being sampled. The default angle is more west than north, so the orange
  # nodes read as external samples rather than decorations on the box corners.
  edges$x_start[sampling] <- edges$x_from[sampling] - edges$width_from[sampling] / 2
  edges$y_start[sampling] <- edges$y_from[sampling] + edges$height_from[sampling] * 0.32
  edges$x_end[sampling] <- edges$x_to[sampling]
  edges$y_end[sampling] <- edges$y_to[sampling]

  edges
}

draw_reproduction_edges <- function(edges) {
  input_edges <- edges[edges$edge_type == "reproduction_input", ]
  output_edges <- edges[edges$edge_type == "birth_output", ]

  list(
    ggplot2::geom_segment(
      data = input_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "#d9368b",
      alpha = 0.28,
      linewidth = 0.95,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
    ),
    ggplot2::geom_segment(
      data = output_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "#d9368b",
      alpha = 0.92,
      linewidth = 2.4,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
    )
  )
}

draw_sampling_edges <- function(edges) {
  sampling_edges <- edges[edges$edge_type == "sampling", ]

  ggplot2::geom_segment(
    data = sampling_edges,
    ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "#e69022",
    linewidth = 0.8,
    lineend = "round",
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  )
}

draw_population_nodes <- function(nodes, overlay) {
  if (nrow(nodes) == 0) {
    return(ggplot2::geom_blank())
  }

  nodes$xmin <- nodes$x - nodes$width / 2
  nodes$xmax <- nodes$x + nodes$width / 2
  nodes$ymin <- nodes$y - nodes$height / 2
  nodes$ymax <- nodes$y + nodes$height / 2
  nodes$fill_value <- overlay_value(nodes, overlay)

  list(
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill_value
      ),
      color = "black",
      linewidth = 0.9
    ),
    ggplot2::scale_fill_gradient(
      low = "#f8f8f4",
      high = "#80b1d3",
      na.value = "#f8f8f4"
    ),
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y + 0.22, label = sex_symbol(sex)),
      size = 8,
      family = "",
      lineheight = 0.9
    ),
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y - 0.30, label = label),
      size = 3.1,
      lineheight = 0.85
    )
  )
}

draw_reproduction_nodes <- function(nodes) {
  if (nrow(nodes) == 0) {
    return(ggplot2::geom_blank())
  }

  nodes$xmin <- nodes$x - nodes$width / 2
  nodes$xmax <- nodes$x + nodes$width / 2
  nodes$ymin <- nodes$y - nodes$height / 2
  nodes$ymax <- nodes$y + nodes$height / 2

  list(
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#222222",
      color = "black",
      linewidth = 0.7
    ),
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label),
      color = "white",
      size = 3.0
    )
  )
}

draw_sample_nodes <- function(nodes) {
  if (nrow(nodes) == 0) {
    return(ggplot2::geom_blank())
  }
  nodes$display_label <- ifelse(nodes$label == "sample", "n", nodes$label)

  list(
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      size = 8,
      stroke = 0.8,
      fill = "#f6a23a",
      color = "black"
    ),
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = display_label),
      size = 3.2,
      fontface = "bold"
    )
  )
}

draw_offspring_pool <- function(shape_data) {
  list(
    ggplot2::geom_polygon(
      data = shape_data$body,
      ggplot2::aes(x = x, y = y, group = node_id),
      fill = "#fff7e8",
      color = "black",
      linewidth = 0.9
    ),
    ggplot2::geom_path(
      data = shape_data$rim,
      ggplot2::aes(x = x, y = y, group = node_id),
      color = "black",
      linewidth = 1.0,
      lineend = "round"
    ),
    ggplot2::geom_text(
      data = shape_data$labels,
      ggplot2::aes(x = x, y = y, label = label),
      size = 3.1
    )
  )
}

make_offspring_shapes <- function(nodes) {
  pools <- nodes[nodes$node_type == "offspring_pool", ]

  if (nrow(pools) == 0) {
    empty <- data.frame(node_id = character(), x = numeric(), y = numeric())
    return(list(body = empty, rim = empty, labels = empty))
  }

  bodies <- list()
  rims <- list()
  labels <- list()

  for (i in seq_len(nrow(pools))) {
    node <- pools[i, ]
    half_width <- node$width / 2
    half_height <- node$height / 2
    left <- node$x - half_width
    right <- node$x + half_width
    top <- node$y + half_height
    bottom <- node$y - half_height

    # A simple pot/crock: flat-ish rim, slightly bowed sides, rounded bottom.
    theta <- seq(205, 335, length.out = 24) * pi / 180
    bottom_arc <- data.frame(
      x = node$x + half_width * 0.86 * cos(theta),
      y = bottom + half_height * 0.72 + half_height * 0.70 * sin(theta)
    )

    body <- rbind(
      data.frame(x = left, y = top),
      data.frame(x = left + 0.08, y = node$y - 0.12),
      bottom_arc,
      data.frame(x = right - 0.08, y = node$y - 0.12),
      data.frame(x = right, y = top)
    )
    body$node_id <- node$node_id

    rim <- data.frame(
      node_id = node$node_id,
      x = seq(left - 0.07, right + 0.07, length.out = 30)
    )
    rim$y <- top + 0.035 * sin(seq(0, pi, length.out = nrow(rim)))

    label <- data.frame(
      node_id = node$node_id,
      x = node$x,
      y = node$y - 0.12,
      label = node$label,
      stringsAsFactors = FALSE
    )

    bodies[[i]] <- body[, c("node_id", "x", "y")]
    rims[[i]] <- rim
    labels[[i]] <- label
  }

  list(
    body = do.call(rbind, bodies),
    rim = do.call(rbind, rims),
    labels = do.call(rbind, labels)
  )
}

make_mini_panels <- function(model, nodes) {
  panels <- list()
  idx <- 1

  for (i in seq_len(nrow(nodes))) {
    node <- nodes[i, ]
    covariates <- model$node_data[[node$node_id]]$covariates

    if (is.null(covariates$length)) {
      next
    }

    curve <- covariates$length
    curve$x_scaled <- scale_to_range(
      curve$x,
      node$x - node$width * 0.34,
      node$x + node$width * 0.34
    )
    curve$y_scaled <- scale_to_range(
      curve$density,
      node$y - node$height * 0.22,
      node$y + node$height * 0.10
    )
    curve$node_id <- node$node_id
    panels[[idx]] <- curve
    idx <- idx + 1
  }

  if (length(panels) == 0) {
    return(data.frame())
  }

  do.call(rbind, panels)
}

draw_mini_panels <- function(panels, overlay) {
  if (overlay != "length_distribution" || nrow(panels) == 0) {
    return(ggplot2::geom_blank())
  }

  ggplot2::geom_path(
    data = panels,
    ggplot2::aes(x = x_scaled, y = y_scaled, group = node_id),
    color = "#4c4c4c",
    linewidth = 0.7
  )
}

overlay_value <- function(nodes, overlay) {
  if (overlay == "abundance") {
    return(nodes$abundance)
  }

  if (overlay == "fecundity") {
    return(nodes$fecundity)
  }

  rep(NA_real_, nrow(nodes))
}

sex_symbol <- function(sex) {
  ifelse(sex == "F", "\u2640", ifelse(sex == "M", "\u2642", ""))
}

scale_to_range <- function(x, lower, upper) {
  rng <- range(x, finite = TRUE)
  if (diff(rng) == 0) {
    return(rep(mean(c(lower, upper)), length(x)))
  }
  lower + (x - rng[1]) / diff(rng) * (upper - lower)
}
