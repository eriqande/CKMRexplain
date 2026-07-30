# Age-time CKMR model plots ----------------------------------------------------
#
# This file is a first stress test for the visual grammar. It is not yet a
# final layout engine. The goal is to see whether the same symbols still read
# when we have multiple time points and multiple non-offspring ages.

plot_ckmr_age_time <- function(model,
                               overlay = c("none", "abundance", "fecundity"),
                               show_reproduction_inputs = TRUE,
                               highlight_edges = attr(model, "highlight_edges", exact = TRUE),
                               path_highlights = attr(model, "path_highlights", exact = TRUE),
                               selected_paths = attr(model, "selected_paths", exact = TRUE)) {
  stopifnot_ckmr_model(model)
  overlay <- match.arg(overlay)

  arranged <- layout_ckmr_age_time(model)
  nodes <- arranged$nodes
  edges <- arranged$edges
  if (is.null(highlight_edges)) {
    highlight_edges <- character()
  }
  edges$highlight <- edges$edge_id %in% highlight_edges
  path_edges <- make_age_time_path_edges(edges, path_highlights, selected_paths)

  population_nodes <- nodes[nodes$node_type == "population", ]
  reproduction_nodes <- nodes[nodes$node_type == "reproduction", ]
  sample_nodes <- nodes[nodes$node_type == "sample", ]

  if (!show_reproduction_inputs) {
    edges <- edges[edges$edge_type != "reproduction_input", ]
  }

  ggplot2::ggplot() +
    draw_age_time_survival_edges(edges) +
    draw_age_time_recruitment_edges(edges) +
    draw_age_time_reproduction_edges(edges) +
    draw_age_time_path_highlights(path_edges) +
    draw_sampling_edges(edges) +
    draw_age_time_population_nodes(population_nodes, overlay) +
    draw_reproduction_nodes(reproduction_nodes) +
    draw_offspring_pool(arranged$offspring_shapes) +
    draw_sample_nodes(sample_nodes) +
    draw_age_time_axes(arranged$axis_labels) +
    ggplot2::coord_equal(
      xlim = arranged$xlim,
      ylim = arranged$ylim,
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.margin = grid::unit(c(8, 14, 8, 28), "pt")
    )
}

layout_ckmr_age_time <- function(model) {
  nodes <- model$nodes
  edges <- model$edges

  population_nodes <- nodes[nodes$node_type == "population", ]
  times <- sort(unique(population_nodes$time))
  ages <- sort(unique(population_nodes$age))

  n_times <- length(times)
  n_ages <- length(ages)

  # The layout gets more compact as the number of time points grows. This is a
  # first approximation; later this should be an explicit layout option.
  dense <- n_times > 5
  col_gap <- if (dense) 1.45 else 2.15
  sex_gap <- if (dense) 0.56 else 0.82
  row_gap <- if (dense) 0.62 else 0.82
  node_size <- if (dense) 0.32 else 0.48
  pot_width <- if (dense) 0.72 else 0.95
  pot_height <- if (dense) 0.48 else 0.62
  repro_width <- if (dense) 0.55 else 0.72
  repro_height <- if (dense) 0.25 else 0.32
  sample_size <- if (dense) 0.20 else 0.26

  time_x <- stats::setNames(seq_along(times) - 1, times) * col_gap
  age_y <- stats::setNames(seq_along(ages) - 1, ages) * row_gap + 1.65

  coords <- data.frame(
    node_id = character(),
    x = numeric(),
    y = numeric(),
    width = numeric(),
    height = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(nodes))) {
    node <- nodes[i, ]

    if (node$node_type == "population") {
      sex_offset <- if (node$sex == "F") -sex_gap / 2 else sex_gap / 2
      coords <- rbind(coords, data.frame(
        node_id = node$node_id,
        x = time_x[as.character(node$time)] + sex_offset,
        y = age_y[as.character(node$age)],
        width = node_size,
        height = node_size,
        stringsAsFactors = FALSE
      ))
    } else if (node$node_type == "reproduction") {
      coords <- rbind(coords, data.frame(
        node_id = node$node_id,
        x = time_x[as.character(node$time)],
        y = 0.95,
        width = repro_width,
        height = repro_height,
        stringsAsFactors = FALSE
      ))
    } else if (node$node_type == "offspring_pool") {
      coords <- rbind(coords, data.frame(
        node_id = node$node_id,
        x = time_x[as.character(node$time)],
        y = 0.08,
        width = pot_width,
        height = pot_height,
        stringsAsFactors = FALSE
      ))
    } else if (node$node_type == "sample") {
      edge <- edges[edges$to == node$node_id & edges$edge_type == "sampling", ][1, ]
      source <- coords[coords$node_id == edge$from, ][1, ]
      coords <- rbind(coords, data.frame(
        node_id = node$node_id,
        x = source$x - node_size * 1.35,
        y = source$y + node_size * 0.82,
        width = sample_size,
        height = sample_size,
        stringsAsFactors = FALSE
      ))
    }
  }

  nodes <- merge(nodes, coords, by = "node_id", all.x = TRUE, sort = FALSE)
  nodes <- nodes[match(model$nodes$node_id, nodes$node_id), ]

  edges <- merge(
    edges,
    nodes[, c("node_id", "x", "y", "width", "height", "node_type", "sex", "age")],
    by.x = "from",
    by.y = "node_id",
    all.x = TRUE,
    sort = FALSE
  )
  names(edges)[names(edges) %in% c("x", "y", "width", "height", "node_type", "sex", "age")] <-
    c("x_from", "y_from", "width_from", "height_from", "type_from", "sex_from", "age_from")

  edges <- merge(
    edges,
    nodes[, c("node_id", "x", "y", "width", "height", "node_type", "sex", "age")],
    by.x = "to",
    by.y = "node_id",
    all.x = TRUE,
    sort = FALSE
  )
  names(edges)[names(edges) %in% c("x", "y", "width", "height", "node_type", "sex", "age")] <-
    c("x_to", "y_to", "width_to", "height_to", "type_to", "sex_to", "age_to")

  edges <- add_age_time_edge_anchors(edges)

  axis_labels <- make_age_time_axis_labels(times, ages, time_x, age_y)

  list(
    nodes = nodes,
    edges = edges,
    offspring_shapes = make_offspring_shapes(nodes),
    axis_labels = axis_labels,
    xlim = range(time_x) + c(-1.0, 1.0),
    ylim = c(-0.42, max(age_y) + 0.70)
  )
}

add_age_time_edge_anchors <- function(edges) {
  edges$x_start <- edges$x_from
  edges$y_start <- edges$y_from
  edges$x_end <- edges$x_to
  edges$y_end <- edges$y_to

  survival <- edges$edge_type == "survival"
  recruitment <- edges$edge_type == "recruitment"
  sampling <- edges$edge_type == "sampling"
  birth_output <- edges$edge_type == "birth_output"
  repro_input <- edges$edge_type == "reproduction_input"

  edges$x_start[survival] <- edges$x_from[survival] + edges$width_from[survival] / 2
  edges$y_start[survival] <- edges$y_from[survival]
  edges$x_end[survival] <- edges$x_to[survival] - edges$width_to[survival] / 2
  edges$y_end[survival] <- edges$y_to[survival]

  # Recruitment carries offspring from the urn at time t into the age-1 male
  # and female boxes at time t + 1. These start from the same point inside the
  # urn so they can later be interpreted as a sex-ratio split.
  edges$x_start[recruitment] <- edges$x_from[recruitment]
  edges$y_start[recruitment] <- edges$y_from[recruitment] + edges$height_from[recruitment] * 0.05
  edges$x_end[recruitment] <- edges$x_to[recruitment]
  edges$y_end[recruitment] <- edges$y_to[recruitment] - edges$height_to[recruitment] / 2

  edges$x_start[repro_input] <- edges$x_from[repro_input]
  edges$y_start[repro_input] <- edges$y_from[repro_input] - edges$height_from[repro_input] / 2
  edges$x_end[repro_input] <- edges$x_to[repro_input]
  edges$y_end[repro_input] <- edges$y_to[repro_input] + edges$height_to[repro_input] / 2

  repro_population <- repro_input & edges$type_from == "population"
  if (any(repro_population)) {
    max_age <- max(edges$age_from[repro_population], na.rm = TRUE)
    age_midpoint <- (min(edges$age_from[repro_population], na.rm = TRUE) + max_age) / 2
    age_offset <- (edges$age_from[repro_population] - age_midpoint) / max(1, max_age)
    sex_direction <- ifelse(edges$sex_from[repro_population] == "F", -1, 1)
    lane_offset <- sex_direction * age_offset * 0.12
    edges$x_start[repro_population] <- edges$x_start[repro_population] + lane_offset
    edges$x_end[repro_population] <- edges$x_end[repro_population] + lane_offset * 0.35
  }

  edges$x_start[birth_output] <- edges$x_from[birth_output]
  edges$y_start[birth_output] <- edges$y_from[birth_output] - edges$height_from[birth_output] / 2
  edges$x_end[birth_output] <- edges$x_to[birth_output]
  edges$y_end[birth_output] <- edges$y_to[birth_output] + edges$height_to[birth_output] / 2

  edges$x_start[sampling] <- edges$x_from[sampling] - edges$width_from[sampling] / 2
  edges$y_start[sampling] <- edges$y_from[sampling] + edges$height_from[sampling] * 0.32
  edges$x_end[sampling] <- edges$x_to[sampling]
  edges$y_end[sampling] <- edges$y_to[sampling]

  edges
}

draw_age_time_population_nodes <- function(nodes, overlay) {
  nodes$xmin <- nodes$x - nodes$width / 2
  nodes$xmax <- nodes$x + nodes$width / 2
  nodes$ymin <- nodes$y - nodes$height / 2
  nodes$ymax <- nodes$y + nodes$height / 2
  nodes$fill_value <- overlay_value(nodes, overlay)
  nodes$sex_label <- nodes$sex

  list(
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill_value
      ),
      color = "black",
      linewidth = 0.45
    ),
    ggplot2::scale_fill_gradient(
      low = "#f8f8f4",
      high = "#80b1d3",
      na.value = "#f8f8f4"
    ),
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = sex_label),
      size = 2.0,
      fontface = "bold"
    )
  )
}

draw_age_time_survival_edges <- function(edges) {
  survival_edges <- edges[edges$edge_type == "survival", ]
  base_edges <- survival_edges[!survival_edges$highlight, ]
  highlight_edges <- survival_edges[survival_edges$highlight, ]

  list(
    ggplot2::geom_segment(
      data = base_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "#2f6df6",
      linewidth = 0.45,
      alpha = 0.70,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
    ),
    ggplot2::geom_segment(
      data = highlight_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "black",
      linewidth = 1.65,
      alpha = 0.95,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.075, "inches"), type = "closed")
    )
  )
}

draw_age_time_recruitment_edges <- function(edges) {
  recruitment_edges <- edges[edges$edge_type == "recruitment", ]
  base_edges <- recruitment_edges[!recruitment_edges$highlight, ]
  highlight_edges <- recruitment_edges[recruitment_edges$highlight, ]

  list(
    ggplot2::geom_segment(
      data = base_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "#7a4cc2",
      linewidth = 0.65,
      alpha = 0.78,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.060, "inches"), type = "closed")
    ),
    ggplot2::geom_segment(
      data = highlight_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "black",
      linewidth = 2.05,
      alpha = 0.96,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.085, "inches"), type = "closed")
    )
  )
}

draw_age_time_reproduction_edges <- function(edges) {
  repro_edges <- edges[edges$edge_type %in% c("reproduction_input", "birth_output"), ]
  repro_edges$edge_alpha <- ifelse(repro_edges$edge_type == "reproduction_input", 0.18, 0.85)
  repro_edges$edge_width <- ifelse(repro_edges$edge_type == "reproduction_input", 0.22, 0.70)

  ggplot2::geom_segment(
    data = repro_edges,
    ggplot2::aes(
      x = x_start, y = y_start, xend = x_end, yend = y_end,
      alpha = edge_alpha, linewidth = edge_width
    ),
    color = "#d9368b",
    lineend = "round",
    arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
  )
}

make_age_time_path_edges <- function(edges, path_highlights, selected_paths) {
  if (is.null(path_highlights) || nrow(path_highlights) == 0) {
    return(data.frame())
  }

  if (is.null(selected_paths)) {
    selected_paths <- unique(path_highlights$path_id)
  }

  path_highlights <- path_highlights[path_highlights$path_id %in% selected_paths, ]
  if (nrow(path_highlights) == 0) {
    return(data.frame())
  }
  if (!"linetype" %in% names(path_highlights)) {
    path_highlights$linetype <- "solid"
  }

  path_highlights$.path_order <- seq_len(nrow(path_highlights))
  path_edges <- merge(
    path_highlights,
    edges,
    by = "edge_id",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  path_edges[order(path_edges$.path_order), ]
}

draw_age_time_path_highlights <- function(edges) {
  if (nrow(edges) == 0) {
    return(ggplot2::geom_blank())
  }

  birth_edges <- edges[edges$edge_role == "birth_output", ]
  parent_edges <- edges[edges$edge_role != "birth_output", ]

  list(
    ggplot2::geom_segment(
      data = birth_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "#b51f73",
      linewidth = 1.85,
      alpha = 0.94,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.070, "inches"), type = "closed")
    ),
    ggplot2::geom_segment(
      data = parent_edges,
      ggplot2::aes(
        x = x_start, y = y_start, xend = x_end, yend = y_end,
        color = color,
        linetype = linetype
      ),
      linewidth = 1.55,
      alpha = 0.96,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.070, "inches"), type = "closed")
    ),
    ggplot2::scale_color_identity(),
    ggplot2::scale_linetype_identity()
  )
}

draw_age_time_axes <- function(axis_labels) {
  list(
    ggplot2::geom_text(
      data = axis_labels$times,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.8
    ),
    ggplot2::geom_text(
      data = axis_labels$ages,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.8,
      hjust = 1
    ),
    ggplot2::geom_text(
      data = axis_labels$sexes,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.4,
      fontface = "bold"
    )
  )
}

make_age_time_axis_labels <- function(times, ages, time_x, age_y) {
  first_x <- min(time_x)
  data_times <- data.frame(
    x = unname(time_x),
    y = -0.30,
    label = paste0("t", times),
    stringsAsFactors = FALSE
  )

  data_ages <- data.frame(
    x = first_x - 0.58,
    y = unname(age_y),
    label = paste0("age ", ages),
    stringsAsFactors = FALSE
  )

  data_sexes <- data.frame(
    x = c(first_x - 0.21, first_x + 0.21),
    y = max(age_y) + 0.42,
    label = c("F", "M"),
    stringsAsFactors = FALSE
  )

  list(times = data_times, ages = data_ages, sexes = data_sexes)
}
