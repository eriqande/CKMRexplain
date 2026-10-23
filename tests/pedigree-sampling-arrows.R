library(CKMRexplain)

sample_times <- tibble::tibble(indiv = c("A", "B-4"), cohort = c(56, 44), age = c(5, 3))
plot <- draw_pedigree_and_cohorts(sample_times = sample_times,
  cohortBarLinewidth = 0.4, cohortBarBlackText = TRUE)
arrow_layers <- function(p) Filter(function(layer)
  inherits(layer$geom, "GeomSegment") && !is.null(layer$geom_params$arrow), p$layers)
layer <- arrow_layers(plot)[[1]]
arrows <- layer$data
layout <- layout_pedigree_and_cohorts()
cells <- layout$cohort_bars
A <- cells[cells$id == "A:1", ]
B <- cells[cells$id == "B:4" & cells$cohort == 44, ]
stopifnot(length(arrow_layers(plot)) == 1L,
  identical(layer$aes_params$colour, "orange"),
  isTRUE(all.equal(layer$aes_params$linewidth, 1.4 * 0.4)),
  identical(arrows$sampling_cohort, c(61, 47)),
  identical(arrows$x, c(A$xmid, B$xmid)),
  identical(arrows$y, layout$BarLate - c(56, 44)),
  identical(arrows$yend, layout$BarLate - c(61, 47)),
  identical(arrows$xend - arrows$x, c(-1.5, 1.5)),
  arrows$xend[1] == A$xlo - (A$xhi - A$xlo),
  arrows$xend[2] == B$xhi + (B$xhi - B$xlo))

custom <- draw_pedigree_and_cohorts(sample_times = sample_times, CohortBoxWidth = 2,
  samplingArrowLinewidth = 0.9, samplingArrowColour = "darkorange")
stopifnot(identical(arrow_layers(custom)[[1]]$data$xend - arrow_layers(custom)[[1]]$data$x,
  c(-3, 3)), identical(arrow_layers(custom)[[1]]$aes_params$linewidth, 0.9))
late <- draw_pedigree_and_cohorts(sample_times = data.frame(indiv = "A-1", cohort = 56, age = 20))
stopifnot(late$coordinates$limits$y[1] <= layout$BarLate - 76 - 0.5,
  any(vapply(late$layers, function(l) "cohort" %in% names(l$data) &&
    76 %in% l$data$cohort, logical(1))))
zero <- draw_pedigree_and_cohorts(sample_times = data.frame(indiv = "B", cohort = 60, age = 0))
stopifnot(arrow_layers(zero)[[1]]$data$y == arrow_layers(zero)[[1]]$data$yend,
  length(arrow_layers(draw_pedigree_and_cohorts())) == 0,
  length(arrow_layers(draw_pedigree_and_cohorts(sample_times = sample_times[FALSE, ]))) == 0)
for (bad in list(list(), data.frame(indiv = "A"),
  data.frame(indiv = "A-0", cohort = 56, age = 5),
  data.frame(indiv = "B-40", cohort = 44, age = 3),
  data.frame(indiv = "A", cohort = 57, age = 5),
  data.frame(indiv = "A", cohort = NA, age = 5),
  data.frame(indiv = "A", cohort = 56, age = -1),
  data.frame(indiv = "A", cohort = 56, age = 1.5))) {
  stopifnot(inherits(try(draw_pedigree_and_cohorts(sample_times = bad), silent = TRUE), "try-error"))
}
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
print(plot)
print(late)
grDevices::dev.off()
