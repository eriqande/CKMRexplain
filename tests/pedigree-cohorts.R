library(CKMRexplain)

layout <- layout_pedigree_and_cohorts(xPedSep = 64, xPosA = 100, xPosB = 400,
  sexA = "M", sexB = "F")
A <- layout$nodes[layout$nodes$tree == "A", ]
stopifnot(nrow(layout$nodes) == 30L, nrow(layout$edges) == 28L,
  identical(A$label, c("A", as.character(2:15))),
  identical(A$x[1:7], c(100, 36, 164, 4, 68, 132, 196)),
  A$cohort_min[4] == 36, A$cohort_max[4] == 50,
  A$cohort_min[7] == 32, A$cohort_max[7] == 44,
  A$y[1] == layout$BarLate - 56,
  all(abs(layout$nodes$ymax[layout$nodes$generation == 3] -
    (layout$BarLate - layout$BarEarly + 0.5)) < 1e-12))

# Enumerate reproductive-age combinations independently of the interval recursion.
for (i in seq_len(nrow(A))) {
  path <- strsplit(A$path[i], "", fixed = TRUE)[[1]]
  ages <- lapply(path, function(sex) if (sex == "M") 3:10 else 6:12)
  possible <- if (length(ages)) {
    sort(unique(56 - rowSums(expand.grid(ages))))
  } else 56
  bar <- layout$cohort_bars[layout$cohort_bars$id == A$id[i], ]
  stopifnot(identical(as.numeric(bar$cohort), as.numeric(possible)),
    all(bar$ymid == layout$BarLate - bar$cohort))
}

for (depth in 0:4) {
  arranged <- layout_pedigree_and_cohorts(n = depth, cA = -4, cB = 2)
  stopifnot(nrow(arranged$nodes) == 2 * (2^(depth+1) - 1),
    nrow(arranged$edges) == nrow(arranged$nodes) - 2,
    !anyDuplicated(arranged$nodes$id))
  for (tree in c("A", "B")) {
    nodes <- arranged$nodes[arranged$nodes$tree == tree, ]
    stopifnot(!anyDuplicated(nodes$x),
      all(diff(sort(nodes$x)) > nodes$xmax[1] - nodes$xmin[1]))
  }
  if (nrow(arranged$edges)) {
    child <- match(arranged$edges$child_id, arranged$nodes$id)
    parent <- match(arranged$edges$parent_id, arranged$nodes$id)
    stopifnot(all(arranged$edges$y == arranged$nodes$ymax[child]),
      all(arranged$edges$yend == arranged$nodes$ymin[parent]),
      all(arranged$edges$yend > arranged$edges$y))
  }
}

bad <- list(list(n = -1), list(n = 1.5), list(Mm = 11), list(cA = NA_real_),
  list(sexB = "unknown"), list(nodeHeight = 0), list(SampleSepFactor = 0),
  list(SampleSepFactor = NA_real_))
for (args in bad) {
  stopifnot(inherits(try(do.call(layout_pedigree_and_cohorts, args), silent = TRUE), "try-error"))
}

# Separation is measured between oldest-ancestor centers, not between roots.
for (depth in 1:4) for (factor in c(1, 2, 3.5)) {
  for (positions in list(list(), list(xPosA = 100), list(xPosB = 400))) {
    arranged <- do.call(layout_pedigree_and_cohorts,
      c(list(n = depth, SampleSepFactor = factor), positions))
    oldest <- arranged$nodes[arranged$nodes$generation == depth, ]
    A <- sort(oldest$x[oldest$tree == "A"])
    B <- sort(oldest$x[oldest$tree == "B"])
    stopifnot(abs((min(B) - max(A)) / diff(A)[1] - factor) < 1e-12)
    if (!is.null(positions$xPosA)) {
      stopifnot(arranged$nodes$x[arranged$nodes$id == "A:1"] == positions$xPosA)
    }
    if (!is.null(positions$xPosB)) {
      stopifnot(arranged$nodes$x[arranged$nodes$id == "B:1"] == positions$xPosB)
    }
  }
}
manual <- layout_pedigree_and_cohorts(xPosA = 100, xPosB = 400, SampleSepFactor = 7)
stopifnot(manual$nodes$x[manual$nodes$id == "A:1"] == 100,
  manual$nodes$x[manual$nodes$id == "B:1"] == 400)

cohort_plot <- draw_pedigree_and_cohorts(showCohorts = TRUE, cohortBarBlackText = FALSE)
text_layers <- Filter(function(layer) inherits(layer$geom, "GeomFitText") &&
  "cohort" %in% names(layer$data) && "generation" %in% names(layer$data), cohort_plot$layers)
colours <- do.call(rbind, lapply(text_layers, function(layer) {
  data.frame(sex = layer$data$sex[1], generation = layer$data$generation[1],
    colour = layer$aes_params$colour)
}))
for (sex in c("M", "F")) {
  shades <- colours[colours$sex == sex & colours$generation > 0, ]
  stopifnot(nrow(shades) == 3, !anyDuplicated(shades$colour),
    all(diff(colMeans(grDevices::col2rgb(shades$colour))) > 0))
}
brewer <- RColorBrewer::brewer.pal(6, "RdBu")
stopifnot(identical(colours$colour[colours$sex == "M" & colours$generation > 0], rev(brewer[4:6])),
  identical(colours$colour[colours$sex == "F" & colours$generation > 0], brewer[1:3]))
bold_plot <- draw_pedigree_and_cohorts(showCohorts = TRUE, cohortBarBold = TRUE,
  cohortBarLinewidth = 0.3)
bar_layers <- Filter(function(layer) all(c("cohort", "xlo") %in% names(layer$data)), bold_plot$layers)
for (layer in bar_layers) {
  if (inherits(layer$geom, "GeomFitText")) {
    stopifnot(layer$aes_params$fontface == "bold")
  } else {
    stopifnot(layer$aes_params$linewidth == 0.3)
  }
}
rectangles <- Filter(function(layer) inherits(layer$geom, "GeomRect"), bold_plot$layers)
stopifnot(all(vapply(tail(rectangles, -1), function(layer) {
  layer$aes_params$linewidth == 0.3
}, logical(1))))
for (args in list(list(generationPalette = "not-a-palette"),
                 list(cohortBarBold = NA), list(cohortBarLinewidth = -1))) {
  stopifnot(inherits(try(do.call(draw_pedigree_and_cohorts, args), silent = TRUE), "try-error"))
}
stopifnot(inherits(draw_pedigree_and_cohorts(showCohorts=TRUE,
  generationPalette=NULL, generationFade=0), "ggplot"))
for (depth in c(1, 4, 6)) {
  deep <- draw_pedigree_and_cohorts(n=depth, showCohorts=TRUE, CohortBoxWidth=0.2,
    cohortBarBlackText=FALSE)
  text <- Filter(function(layer) inherits(layer$geom, "GeomFitText") &&
    "generation" %in% names(layer$data) && "cohort" %in% names(layer$data), deep$layers)
  stopifnot(length(unique(vapply(text,function(layer) layer$aes_params$colour,character(1)))) == 2*depth)
}
plain <- draw_pedigree_and_cohorts()
stopifnot(inherits(plain$layers[[1]]$geom, "GeomSegment"),
  nrow(plain$layers[[1]]$data) == 30,
  plain$layers[[1]]$aes_params$linetype == "dashed",
  plain$layers[[1]]$aes_params$colour == "black",
  cohort_plot$layers[[1]]$aes_params$colour == "grey70")
no_guides <- draw_layer_pedigree(layout_pedigree_and_cohorts(), showGuides = FALSE)
stopifnot(is.null(no_guides[[1]]))

# Render both views and the zero-generation, all-male edge case to a vector device.
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
for (plot in list(draw_pedigree_and_cohorts(),
                 draw_pedigree_and_cohorts(showCohorts = TRUE),
                 bold_plot,
                 draw_pedigree_and_cohorts(n = 0, sexA = "M", sexB = "M", showGrid = FALSE))) {
  stopifnot(inherits(plot, "ggplot"))
  print(plot)
}
grDevices::dev.off()
