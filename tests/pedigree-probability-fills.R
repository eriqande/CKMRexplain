library(CKMRexplain)

ages <- 1:3
spd <- list("max-age" = 3, "male-prob-repro" = c(1, 0, 1),
  "fem-prob-repro" = c(0, 1, 1), "male-asrp" = ages^2,
  "fem-asrf" = ages^2, "male-surv-probs" = rep(0.85, 3),
  "fem-surv-probs" = rep(0.85, 3))
census <- expand.grid(year = 0:20, age = ages, pop = 0)
census$male <- 10 + census$year + census$age
census$female <- 20 + 2 * census$year + census$age
demography <- cohort_demography_from_spd(spd, census)
probability_plot <- draw_pedigree_and_cohorts(demography = demography,
  cA = 15, cB = 17, n = 2, Am = 99, Af = 99, Mm = 99, Mf = 99)
fill_layers <- function(p) Filter(function(l)
  inherits(l$geom, "GeomRect") && !is.null(l$mapping$fill), p$layers)
cells <- do.call(rbind, lapply(fill_layers(probability_plot), function(l) l$data))
expected <- pedigree_cohort_probabilities(demography, 15, 17, 2)$cohort_probabilities
keys <- function(x) paste(x$tree, x$index, x$cohort)
stopifnot(nrow(cells) == nrow(expected), !anyDuplicated(keys(cells)),
  identical(cells$probability, expected$probability[match(keys(cells), keys(expected))]),
  any(cells$probability == 0), any(cells$probability == 1),
  all(cells$cohort[cells$id == "A:2"] == 12:14))
plain <- draw_pedigree_and_cohorts(showCohorts = TRUE, cA = 15, cB = 17,
  n = 2, Am = 3, Af = 3, Mm = 1, Mf = 2)
segments <- function(p) unname(lapply(Filter(function(l)
  inherits(l$geom, "GeomSegment"), p$layers), function(l) l$aes_params$colour))
stopifnot(identical(segments(plain), segments(probability_plot)),
  length(fill_layers(plain)) == 0,
  probability_plot$theme$legend.position == "bottom",
  plain$theme$legend.position == "none",
  !any(vapply(fill_layers(probability_plot), function(l) any(is.na(l$data$index)), logical(1))))
focused <- draw_pedigree_and_cohorts(demography = demography, cA = 15, cB = 17, n = 2,
  focal_indiv_A = "A-2", focal_cohort_A = 12)
static_fills <- Filter(function(l) inherits(l$geom, "GeomRect") &&
  identical(l$aes_params$fill, "#FF00FF"), focused$layers)
stopifnot(length(static_fills) == 1, static_fills[[1]]$data$cohort == 12,
  !any(vapply(focused$layers, function(l) identical(l$aes_params$fill, "grey88"), logical(1))))
built <- ggplot2::ggplot_build(probability_plot)
scale <- built$plot$scales$get_scales("fill")
stopifnot(identical(scale$limits, c(0, 1)),
  identical(scale$trans$name, "sqrt"))
linear <- draw_pedigree_and_cohorts(demography = demography, cA = 15, cB = 17,
  n = 2, probabilityFillTransform = "identity")
linear_scale <- ggplot2::ggplot_build(linear)$plot$scales$get_scales("fill")
stopifnot(identical(linear_scale$trans$name, "identity"),
  identical(linear_scale$limits, c(0, 1)),
  identical(linear_scale$breaks, seq(0, 1, by = 0.2)),
  identical(lapply(fill_layers(linear), function(l) l$data$probability),
    lapply(fill_layers(probability_plot), function(l) l$data$probability)),
  !identical(linear_scale$map(linear_scale$transform(0.1)), scale$map(scale$transform(0.1))))
zero <- draw_pedigree_and_cohorts(demography = demography, cA = 0, cB = 20, n = 0)
stopifnot(all(do.call(rbind, lapply(fill_layers(zero), function(l) l$data))$probability == 1))
hidden <- draw_pedigree_and_cohorts(demography = demography,
  cA = 15, cB = 17, n = 2, showCohorts = FALSE)
stopifnot(length(fill_layers(hidden)) == 0, hidden$theme$legend.position == "none")
grDevices::pdf(tempfile(fileext = ".pdf"), width = 14, height = 9)
print(probability_plot)
print(linear)
print(focused)
print(zero)
grDevices::dev.off()
