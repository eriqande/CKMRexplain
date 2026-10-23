library(CKMRexplain)

spec <- tibble::tibble(focalA = "A", focalB = "B", cohortA = 56, cohortB = 60,
  ancestorA = "A-2", ancestorB = "B-2")
layout <- layout_pedigree_and_cohorts()
compute <- CKMRexplain:::pedigree_shared_H11
auto <- compute(layout, spec)[[1]]
explicit <- spec
explicit$ancestor_cohorts <- list(50:53)
stopifnot(identical(auto, compute(layout, explicit)[[1]]),
  identical(as.numeric(auto$cohorts), as.numeric(50:53)),
  nrow(auto$descent) == 8, nrow(auto$matches) == 4,
  identical(auto$matches$colour, grDevices::rainbow(4)),
  all(auto$matches$y == auto$matches$yend),
  isTRUE(all.equal(auto$survival$y, layout$BarLate - 56 + seq(1/6, -1/6, length.out = 4))),
  isTRUE(all.equal(auto$survival$yend, layout$BarLate - 60 + seq(1/6, -1/6, length.out = 4))),
  identical(auto$survival$colour, auto$matches$colour),
  identical(auto$descent$yend[1:4], auto$survival$y),
  identical(auto$descent$yend[5:8], auto$survival$yend),
  isTRUE(all.equal(auto$rectangles$ymax - auto$rectangles$ymin, rep(1/3, 2))),
  all(auto$rectangles$xmax - auto$rectangles$xmin == 1/4),
  auto$points$x[1] == max(layout$cohort_bars$xhi[layout$cohort_bars$tree == "A"]) + 1,
  auto$points$x[2] == min(layout$cohort_bars$xlo[layout$cohort_bars$tree == "B"]) - 1)
explicit$ancestor_cohorts <- list(NULL)
stopifnot(identical(auto, compute(layout, explicit)[[1]]))
explicit$ancestor_cohorts <- list(c(53, 50))
stopifnot(identical(compute(layout, explicit)[[1]]$cohorts, c(50, 53)))
explicit$ancestor_cohorts <- list(51)
single <- compute(layout, explicit)[[1]]
stopifnot(single$survival$y == layout$BarLate - 56,
  single$survival$yend == layout$BarLate - 60)
older <- spec
older$focalA <- "A-2"; older$focalB <- "B-2"
older$cohortA <- 51; older$cohortB <- 54
older$ancestorA <- "A-4"; older$ancestorB <- "B-4"
stopifnot(identical(as.numeric(compute(layout, older)[[1]]$cohorts), as.numeric(44:48)))
female <- spec
female$ancestorA <- "A-3"; female$ancestorB <- "B-3"
stopifnot(identical(as.numeric(compute(layout, female)[[1]]$cohorts), as.numeric(48:50)),
  length(compute(layout, rbind(spec, female))) == 2,
  length(compute(layout, spec[FALSE, ])) == 0)
bad <- list(spec[, -1], transform(spec, ancestorA = "A-4"),
  transform(spec, ancestorB = "B-3"), transform(spec, focalA = "B"),
  transform(spec, cohortA = 57), transform(spec, cohortA = NA_real_))
for (values in list(integer(), c(50, 50), 49, 50.5, NA_real_, "50")) {
  invalid <- spec
  invalid$ancestor_cohorts <- list(values)
  bad[[length(bad) + 1]] <- invalid
}
for (invalid in bad) stopifnot(inherits(try(compute(layout, invalid), silent = TRUE), "try-error"))
far <- spec; far$cohortB <- 100
stopifnot(inherits(try(compute(layout_pedigree_and_cohorts(cB = 100), far), silent = TRUE), "try-error"))
plot <- draw_pedigree_and_cohorts(shared_H11 = spec, cohortBarBlackText = TRUE,
  focal_indiv_A = "A", focal_cohort_A = 56, focal_indiv_B = "B", focal_cohort_B = 60)
invisible(ggplot2::ggplot_build(plot))
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
print(plot)
print(draw_pedigree_and_cohorts(shared_H11 = older, cohortBarBlackText = TRUE))
print(draw_pedigree_and_cohorts(shared_H11 = female))
grDevices::dev.off()
