library(CKMRexplain)

base <- layout_pedigree_and_cohorts()
blank <- highlight_pedigree_cohorts(base)
stopifnot(!any(blank$cohort_bars$reachable), !any(blank$cohort_bars$focal),
  all(blank$cohort_bars$highlight == "none"), nrow(blank$focal_selections) == 0)

focused <- highlight_pedigree_cohorts(base,
  focal_indiv_A = "A-2", focal_cohort_A = 51,
  focal_indiv_B = "B-6", focal_cohort_B = 44)
bars <- focused$cohort_bars
stopifnot(sum(bars$focal) == 2,
  all(bars$reachable[bars$focal]),
  setequal(unique(bars$id[bars$reachable & bars$tree == "A"]),
    paste0("A:", c(2,4,5,8,9,10,11))),
  setequal(unique(bars$id[bars$reachable & bars$tree == "B"]),
    paste0("B:", c(6,12,13))),
  identical(base$nodes, focused$nodes), identical(base$edges, focused$edges))

# Enumerate the remaining sex paths from each selected cohort independently.
for (tree in c("A", "B")) {
  selection <- focused$focal_selections[focused$focal_selections$tree == tree, ]
  focal_path <- base$nodes$path[base$nodes$id == selection$id]
  for (i in which(base$nodes$tree == tree)) {
    node <- base$nodes[i, ]
    possible <- numeric()
    if (startsWith(node$path, focal_path)) {
      suffix <- substring(node$path, nchar(focal_path) + 1)
      sexes <- strsplit(suffix, "", fixed = TRUE)[[1]]
      ages <- lapply(sexes, function(sex) if (sex == "M") 3:10 else 6:12)
      possible <- if (length(ages)) sort(unique(selection$cohort - rowSums(expand.grid(ages)))) else selection$cohort
    }
    observed <- bars$cohort[bars$id == node$id & bars$reachable]
    stopifnot(identical(as.numeric(observed), as.numeric(possible)))
  }
}

roots <- highlight_pedigree_cohorts(base,
  focal_indiv_A = "A", focal_cohort_A = 56,
  focal_indiv_B = "B", focal_cohort_B = 60)
stopifnot(all(roots$cohort_bars$reachable), sum(roots$cohort_bars$focal) == 2)
alias <- highlight_pedigree_cohorts(base,
  focal_indiv_A = "A-1", focal_cohort_A = 56,
  focal_indiv_B = "B-1", focal_cohort_B = 60)
stopifnot(identical(roots$cohort_bars, alias$cohort_bars))
one <- highlight_pedigree_cohorts(base, focal_indiv_A = "A-15", focal_cohort_A = 30)
stopifnot(sum(one$cohort_bars$reachable) == 1, !any(one$cohort_bars$reachable[one$cohort_bars$tree == "B"]))
reset <- highlight_pedigree_cohorts(focused)
stopifnot(identical(reset$cohort_bars, blank$cohort_bars))
zero <- highlight_pedigree_cohorts(layout_pedigree_and_cohorts(n = 0),
  focal_indiv_A = "A", focal_cohort_A = 56)
stopifnot(sum(zero$cohort_bars$reachable) == 1, sum(zero$cohort_bars$focal) == 1)

bad <- list(list(focal_indiv_A = "A"), list(focal_cohort_A = 56),
  list(focal_indiv_A = "B-2", focal_cohort_A = 51),
  list(focal_indiv_A = "A-0", focal_cohort_A = 51),
  list(focal_indiv_A = "A-16", focal_cohort_A = 30),
  list(focal_indiv_A = "A", focal_cohort_A = 51),
  list(focal_indiv_A = "A-2", focal_cohort_A = 54),
  list(focal_indiv_A = "A-2", focal_cohort_A = 51.5),
  list(focal_indiv_A = "A-2", focal_cohort_A = NA_real_))
for (args in bad) {
  stopifnot(inherits(try(do.call(highlight_pedigree_cohorts, c(list(layout = base), args)),
    silent = TRUE), "try-error"))
}

plot <- draw_pedigree_and_cohorts(focal_indiv_A = "A-2", focal_cohort_A = 51,
  focal_indiv_B = "B-6", focal_cohort_B = 44, cohortBarBold = TRUE,
  cohortBarBlackText = FALSE)
fills <- Filter(function(layer) inherits(layer$geom, "GeomRect") &&
  identical(layer$aes_params$fill, "#FF00FF"), plot$layers)
cells <- do.call(rbind, lapply(fills, function(layer) layer$data))
stopifnot(setequal(paste(cells$id, cells$cohort), c("A:2 51", "B:6 44")))
grey <- Filter(function(layer) inherits(layer$geom, "GeomRect") &&
  identical(layer$aes_params$fill, "grey88"), plot$layers)
cells <- do.call(rbind, lapply(grey, function(layer) layer$data))
expected <- bars[bars$reachable & !bars$focal, ]
stopifnot(setequal(paste(cells$id, cells$cohort), paste(expected$id, expected$cohort)))
labels <- Filter(function(layer) inherits(layer$geom, "GeomFitText") &&
  "focal" %in% names(layer$data) && any(layer$data$focal), plot$layers)
stopifnot(sum(vapply(labels,function(layer) nrow(layer$data),integer(1))) == 2,
  all(vapply(labels,function(layer) layer$aes_params$colour == "black",logical(1))))
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
print(plot)
print(draw_pedigree_and_cohorts(focal_indiv_A="A",focal_cohort_A=56,
  focal_indiv_B="B",focal_cohort_B=60))
grDevices::dev.off()

black_text <- draw_pedigree_and_cohorts(focal_indiv_A = "A-2", focal_cohort_A = 51,
  focal_indiv_B = "B-6", focal_cohort_B = 44, cohortBarBold = TRUE,
  cohortBarBlackText = TRUE)
cohort_labels <- Filter(function(layer) inherits(layer$geom, "GeomFitText") &&
  "focal" %in% names(layer$data), black_text$layers)
stopifnot(length(cohort_labels) > 0,
  all(vapply(cohort_labels, function(layer) layer$aes_params$colour == "black", logical(1))),
  all(vapply(cohort_labels, function(layer) layer$aes_params$fontface == "bold", logical(1))))
outlines <- function(p) lapply(Filter(function(layer)
  inherits(layer$geom, "GeomSegment"), p$layers), function(layer) layer$aes_params$colour)
stopifnot(identical(unname(outlines(plot)), unname(outlines(black_text))))
for (invalid in list(NA, 1, "black", logical(), c(TRUE, FALSE))) {
  stopifnot(inherits(try(draw_pedigree_and_cohorts(cohortBarBlackText = invalid),
    silent = TRUE), "try-error"))
}
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
print(black_text)
grDevices::dev.off()
