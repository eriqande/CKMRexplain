pedigree_sampling_arrows <- function(layout, sample_times) {
  arrows <- data.frame(x = numeric(), y = numeric(), xend = numeric(),
    yend = numeric(), sampling_cohort = numeric())
  if (is.null(sample_times)) return(arrows)
  if (!is.data.frame(sample_times) ||
      !all(c("indiv", "cohort", "age") %in% names(sample_times))) {
    stop("`sample_times` must be a data frame with indiv, cohort, and age columns.", call. = FALSE)
  }
  if (!is.character(sample_times$indiv) || anyNA(sample_times$indiv) ||
      any(!grepl("^[AB](-[1-9][0-9]*)?$", sample_times$indiv))) {
    stop("Use 'A', 'B', or numbered slots such as 'B-4' in `sample_times$indiv`.", call. = FALSE)
  }
  for (column in c("cohort", "age")) {
    values <- sample_times[[column]]
    if (!is.numeric(values) || any(!is.finite(values)) || any(values != floor(values))) {
      stop("`sample_times$", column, "` must contain finite integers.", call. = FALSE)
    }
  }
  if (any(sample_times$age < 0)) {
    stop("`sample_times$age` must be nonnegative.", call. = FALSE)
  }
  bars <- layout$cohort_bars
  for (i in seq_len(nrow(sample_times))) {
    individual <- sample_times$indiv[i]
    tree <- substr(individual, 1, 1)
    index <- if (individual == tree) 1 else as.numeric(sub("^[AB]-", "", individual))
    cell <- bars[bars$tree == tree & bars$index == index &
      bars$cohort == sample_times$cohort[i], ]
    if (nrow(cell) != 1L) {
      stop("No cohort cell for ", individual, " at cohort ", sample_times$cohort[i],
        "; choose an individual and cohort present in the pedigree.", call. = FALSE)
    }
    sampling_cohort <- sample_times$cohort[i] + sample_times$age[i]
    if (!is.finite(sampling_cohort)) stop("Sampling times must be finite.", call. = FALSE)
    arrows[i, ] <- list(cell$xmid, cell$ymid,
      if (tree == "A") cell$xlo - (cell$xhi - cell$xlo) else cell$xhi + (cell$xhi - cell$xlo),
      layout$BarLate - sampling_cohort, sampling_cohort)
  }
  arrows
}
