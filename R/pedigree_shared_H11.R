pedigree_shared_H11 <- function(layout, shared_H11) {
  result <- list()
  if (is.null(shared_H11)) return(result)
  required <- c("focalA", "focalB", "cohortA", "cohortB", "ancestorA", "ancestorB")
  if (!is.data.frame(shared_H11) || !all(required %in% names(shared_H11))) {
    stop("`shared_H11` must be a data frame with focalA, focalB, cohortA, cohortB, ancestorA, and ancestorB.", call. = FALSE)
  }
  if ("ancestor_cohorts" %in% names(shared_H11) && !is.list(shared_H11$ancestor_cohorts)) {
    stop("`ancestor_cohorts` must be a list column (NULL entries select cohorts automatically).", call. = FALSE)
  }
  nodes <- layout$nodes
  bars <- layout$cohort_bars
  resolve <- function(name, tree) {
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !grepl(paste0("^", tree, "(-[1-9][0-9]*)?$"), name)) {
      stop("Use '", tree, "' or '", tree, "-<index>' in the corresponding shared_H11 column.", call. = FALSE)
    }
    index <- if (name == tree) 1 else as.numeric(sub("^[AB]-", "", name))
    node <- nodes[nodes$tree == tree & nodes$index == index, ]
    if (nrow(node) != 1L) stop("Individual ", name, " is not present in the pedigree.", call. = FALSE)
    node
  }
  for (row in seq_len(nrow(shared_H11))) {
    focal <- ancestor <- list()
    cohorts <- numeric(2)
    for (j in seq_along(c("A", "B"))) {
      tree <- c("A", "B")[j]
      focal[[j]] <- resolve(shared_H11[[paste0("focal", tree)]][row], tree)
      ancestor[[j]] <- resolve(shared_H11[[paste0("ancestor", tree)]][row], tree)
      cohort <- shared_H11[[paste0("cohort", tree)]][row]
      if (!is.numeric(cohort) || length(cohort) != 1L || !is.finite(cohort) ||
          cohort != floor(cohort) || cohort < focal[[j]]$cohort_min || cohort > focal[[j]]$cohort_max) {
        stop("Each focal cohort must be an integer in its individual's cohort bar.", call. = FALSE)
      }
      cohorts[j] <- cohort
      if (is.na(ancestor[[j]]$parent_index) || ancestor[[j]]$parent_index != focal[[j]]$index) {
        stop("Each shared ancestor must be an immediate parent of its focal individual.", call. = FALSE)
      }
    }
    if (ancestor[[1]]$sex != ancestor[[2]]$sex) {
      stop("Shared ancestors must have the same sex.", call. = FALSE)
    }
    ages <- layout$reproductive_ages[layout$reproductive_ages$sex == ancestor[[1]]$sex, ]
    low <- max(cohorts - ages$max_age, ancestor[[1]]$cohort_min, ancestor[[2]]$cohort_min)
    high <- min(cohorts - ages$min_age, ancestor[[1]]$cohort_max, ancestor[[2]]$cohort_max)
    if (low > high) stop("No mutually feasible shared ancestor cohorts in row ", row, ".", call. = FALSE)
    possible <- seq.int(low, high)
    selected <- if ("ancestor_cohorts" %in% names(shared_H11)) shared_H11$ancestor_cohorts[[row]] else NULL
    if (is.null(selected)) selected <- possible
    if (!is.numeric(selected) || !length(selected) || any(!is.finite(selected)) ||
        any(selected != floor(selected)) || any(!selected %in% possible) || anyDuplicated(selected)) {
      stop("`ancestor_cohorts` must contain distinct feasible integer cohorts in row ", row, ".", call. = FALSE)
    }
    selected <- sort(selected)
    colours <- grDevices::rainbow(length(selected))
    width <- bars$xhi[1] - bars$xlo[1]
    ends <- c(max(bars$xhi[bars$tree == "A"]) + width,
      min(bars$xlo[bars$tree == "B"]) - width)
    y <- layout$BarLate - cohorts
    offsets <- if (length(selected) == 1L) 0 else seq(1 / 6, -1 / 6, length.out = length(selected))
    descent <- do.call(rbind, lapply(1:2, function(j) data.frame(
      x = ancestor[[j]]$x, y = layout$BarLate - selected,
      xend = ends[j], yend = y[j] + offsets, cohort = selected, colour = colours)))
    matches <- data.frame(x = ancestor[[1]]$x, xend = ancestor[[2]]$x,
      y = layout$BarLate - selected, yend = layout$BarLate - selected,
      cohort = selected, colour = colours)
    focal_links <- data.frame(x = c(focal[[1]]$x, focal[[2]]$x), y = y,
      xend = ends, yend = y)
    survival <- data.frame(x = ends[1], y = y[1] + offsets,
      xend = ends[2], yend = y[2] + offsets, cohort = selected, colour = colours)
    points <- data.frame(x = ends, y = y)
    rectangles <- data.frame(xmin = ends - width / 8, xmax = ends + width / 8,
      ymin = y - 1 / 6, ymax = y + 1 / 6)
    result[[row]] <- list(cohorts = selected, descent = descent, matches = matches,
      focal_links = focal_links, survival = survival, points = points, rectangles = rectangles)
  }
  result
}

pedigree_shared_H11_layers <- function(layout, relationships, linewidth) {
  segments <- function(data, colour, linetype = "solid") {
    ggplot2::geom_segment(data = data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = colour, linetype = linetype, linewidth = linewidth, inherit.aes = FALSE)
  }
  lapply(relationships, function(relation) {
    coloured <- lapply(unique(relation$matches$colour), function(colour) list(
      segments(relation$matches[relation$matches$colour == colour, ], colour, "dotted"),
      segments(relation$descent[relation$descent$colour == colour, ], colour),
      segments(relation$survival[relation$survival$colour == colour, ], colour)))
    list(segments(relation$focal_links, "black", "dotted"), coloured,
      ggplot2::geom_rect(data = relation$rectangles,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "black", colour = NA, inherit.aes = FALSE))
  })
}
