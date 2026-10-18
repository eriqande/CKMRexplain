#' Calculate pedigree positions and feasible ancestral cohorts
#'
#' Index 1 is the sample, 2i is the father of i, and 2i+1 is its mother.
#' Generation heights are diagrammatic; cohort intervals are actual birth years.
#' Every integer reproductive age between the stated limits is considered possible.
#' @param Am,Af Maximum reproductive ages of males and females.
#' @param Mm,Mf Minimum reproductive ages of males and females.
#' @param cA,cB Birth cohorts of samples A and B.
#' @param aA,aB Ages when samples A and B are sampled.
#' @param n Number of ancestral generations; zero draws the samples only.
#' @param xPosBar Horizontal position of the cohort reference bar.
#' @param xPosA,xPosB Sample positions. NULL positions a tree relative to the
#' other tree; if both are NULL the trees are placed to the right of the reference.
#' @param xPedSep Sample-to-parent horizontal distance, halved each generation.
#' @param CohortBoxWidth Width of cohort bars in x units.
#' @param sexA,sexB Sample sexes, either "M" or "F".
#' @param nodeWidth,nodeHeight Node dimensions in plot units. NULL derives them
#' from the horizontal and vertical spacing.
#' @param SampleSepFactor Positive multiple of adjacent oldest-ancestor spacing
#' used for the center-to-center gap between the trees' nearest oldest ancestors.
#' Only used when at least one sample position is NULL. For n = 0, node width
#' replaces oldest-ancestor spacing.
#' @return A list with `nodes`, `edges`, `cohort_bars`, `reference_bar`, `grid`,
#' `reproductive_ages`, `BarEarly`, `BarLate`, `aspect_ratio`, `xlim`, and `ylim`.
#' @md
#' @export
layout_pedigree_and_cohorts <- function(
  Am = 10, Af = 12, Mm = 3, Mf = 6,
  cA = 56, cB = 60, aA = 5, aB = 3, n = 3,
  xPosBar = 0, xPosA = NULL, xPosB = NULL, xPedSep = 16,
  CohortBoxWidth = 1, sexA = "F", sexB = "F",
  nodeWidth = NULL, nodeHeight = NULL, SampleSepFactor = 2
) {
  values <- list(Am = Am, Af = Af, Mm = Mm, Mf = Mf, cA = cA, cB = cB,
    aA = aA, aB = aB, n = n, xPosBar = xPosBar, xPedSep = xPedSep,
    CohortBoxWidth = CohortBoxWidth, SampleSepFactor = SampleSepFactor)
  for (nm in names(values)) {
    x <- values[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop("`", nm, "` must be one finite number.", call. = FALSE)
    }
  }
  integers <- unlist(values[c("Am", "Af", "Mm", "Mf", "cA", "cB", "aA", "aB", "n")])
  if (any(integers != floor(integers)) || n < 0 || min(Mm, Mf) < 1 ||
      Am < Mm || Af < Mf || min(aA, aB) < 0 || xPedSep <= 0 ||
      CohortBoxWidth <= 0 || SampleSepFactor <= 0) {
    stop("Use integer cohorts, ages and n, valid positive reproductive ages, and positive widths.", call. = FALSE)
  }
  if (length(sexA) != 1L || length(sexB) != 1L ||
      !sexA %in% c("M", "F") || !sexB %in% c("M", "F")) {
    stop("`sexA` and `sexB` must each be 'M' or 'F'.", call. = FALSE)
  }
  for (nm in c("xPosA", "xPosB", "nodeWidth", "nodeHeight")) {
    x <- get(nm)
    if (!is.null(x) && (!is.numeric(x) || length(x) != 1L || !is.finite(x))) {
      stop("`", nm, "` must be NULL or one finite number.", call. = FALSE)
    }
  }
  BarEarly <- min(cA, cB) - n * max(Am, Af)
  BarLate <- max(cA + aA, cB + aB)
  top <- BarLate - BarEarly + 0.5
  spacing <- xPedSep / 2^max(0, n - 1)
  if (is.null(nodeWidth)) nodeWidth <- 0.6 * spacing
  if (is.null(nodeHeight)) {
    gap <- (min(cA, cB) - BarEarly) / max(1, n)
    nodeHeight <- if (n == 0) 1.2 else min(1.2, gap / 3)
  }
  if (nodeWidth <= 0 || nodeHeight <= 0) stop("Node dimensions must be positive.", call. = FALSE)
  reach <- 2 * xPedSep * (1 - 2^-n)
  margin <- max(nodeWidth, CohortBoxWidth)
  oldest_spacing <- if (n == 0) nodeWidth else 2 * spacing
  sample_separation <- 2 * reach + SampleSepFactor * oldest_spacing
  if (is.null(xPosA) && is.null(xPosB)) {
    xPosA <- xPosBar + reach + 2 * margin
    xPosB <- xPosA + sample_separation
  } else if (is.null(xPosA)) {
    xPosA <- xPosB - sample_separation
  } else if (is.null(xPosB)) {
    xPosB <- xPosA + sample_separation
  }
  if (n > 0 && max(nodeWidth, CohortBoxWidth) >= spacing) {
    warning("Some nodes or cohort bars may overlap: reduce widths or increase xPedSep.", call. = FALSE)
  }
  make_nodes <- function(tree, cohort, xpos, sample_sex) {
    index <- seq_len(2^(n + 1) - 1)
    generation <- floor(log2(index))
    nodes <- data.frame(tree = tree, index = index, id = paste(tree, index, sep = ":"),
      generation = generation, parent_index = c(NA_integer_, floor(index[-1] / 2)),
      label = c(tree, as.character(index[-1])),
      sex = c(sample_sex, ifelse(index[-1] %% 2 == 0, "M", "F")),
      path = "", x = xpos, y = BarLate - cohort,
      cohort_min = cohort, cohort_max = cohort, stringsAsFactors = FALSE)
    if (n > 0) {
      nodes$y <- (BarLate - cohort) + generation / n *
        (top - nodeHeight / 2 - (BarLate - cohort))
      for (i in index[-1]) {
        child <- nodes$parent_index[i]
        female <- nodes$sex[i] == "F"
        step <- xPedSep / 2^(generation[i] - 1)
        nodes$x[i] <- nodes$x[child] + if (female) step else -step
        nodes$path[i] <- paste0(nodes$path[child], nodes$sex[i])
        nodes$cohort_min[i] <- nodes$cohort_min[child] - if (female) Af else Am
        nodes$cohort_max[i] <- nodes$cohort_max[child] - if (female) Mf else Mm
      }
    }
    nodes$xmin <- nodes$x - nodeWidth / 2
    nodes$xmax <- nodes$x + nodeWidth / 2
    nodes$ymin <- nodes$y - nodeHeight / 2
    nodes$ymax <- nodes$y + nodeHeight / 2
    nodes
  }
  nodes <- rbind(make_nodes("A", cA, xPosA, sexA), make_nodes("B", cB, xPosB, sexB))
  ancestors <- nodes[nodes$generation > 0, ]
  child <- match(paste(ancestors$tree, ancestors$parent_index, sep = ":"), nodes$id)
  edges <- data.frame(child_id = nodes$id[child], parent_id = ancestors$id,
    x = nodes$x[child], y = nodes$ymax[child], xend = ancestors$x, yend = ancestors$ymin)
  make_bar <- function(cohorts, xpos) {
    data.frame(cohort = cohorts, ymid = BarLate - cohorts, xmid = xpos,
      xlo = xpos - CohortBoxWidth / 2, xhi = xpos + CohortBoxWidth / 2)
  }
  cohort_bars <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
    bar <- make_bar(seq.int(nodes$cohort_min[i], nodes$cohort_max[i]), nodes$x[i])
    bar$id <- nodes$id[i]
    bar$tree <- nodes$tree[i]
    bar$index <- nodes$index[i]
    bar$sex <- nodes$sex[i]
    bar$generation <- nodes$generation[i]
    bar
  }))
  reference_bar <- make_bar(seq.int(BarEarly, BarLate), xPosBar)
  grid <- data.frame(cohort = seq.int(BarEarly, BarLate))
  grid$y <- BarLate - grid$cohort
  grid$major <- grid$cohort %% 10 == 0
  grid <- grid[grid$cohort %% 5 == 0, ]
  limits_x <- range(c(nodes$xmin, nodes$xmax, reference_bar$xlo,
    reference_bar$xhi, cohort_bars$xlo, cohort_bars$xhi)) + c(-margin, margin)
  limits_y <- range(c(-0.5, top, nodes$ymin, nodes$ymax)) + c(-0.5, 0.5)
  list(nodes = nodes, edges = edges, cohort_bars = cohort_bars,
    reproductive_ages = data.frame(sex = c("M", "F"),
      min_age = c(Mm, Mf), max_age = c(Am, Af)),
    reference_bar = reference_bar, grid = grid, BarEarly = BarEarly, BarLate = BarLate,
    aspect_ratio = nodeWidth / nodeHeight, xlim = limits_x, ylim = limits_y)
}
