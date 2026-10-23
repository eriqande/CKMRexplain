#' Calculate unconditional birth-cohort distributions for pedigree positions
#'
#' Recursively propagates each child's cohort distribution through sex-specific
#' parental age distributions. Each pedigree position has total probability one;
#' no sharing/non-sharing conditions, survival factors, or focal selections are used.
#' A conservative history check requires the full potential ancestral birth-cohort
#' range to lie within the recorded census window. It does not truncate or
#' renormalize ancestry at the historical boundary.
#' @param demography Inputs returned by cohort_demography_from_spd().
#' @param cA,cB Known birth cohorts of the two samples.
#' @param n Number of ancestral generations, zero or greater.
#' @return A list containing nodes, a long cohort_probabilities table (tree,
#' individual, index, generation, sex, cohort, probability), and demography.
#' Root sex is NA because it does not enter this calculation. Zero-probability
#' cohorts inside the structural minimum--maximum interval are retained.
#' @export
pedigree_cohort_probabilities <- function(demography, cA = 56, cB = 60, n = 3) {
  values <- list(cA = cA, cB = cB, n = n)
  for (name in names(values)) {
    value <- values[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != floor(value)) {
      stop("`", name, "` must be one finite integer.", call. = FALSE)
    }
  }
  if (n < 0 || n > 20) stop("`n` must be between 0 and 20.", call. = FALSE)
  if (!is.list(demography) || !all(c("E", "N", "years", "reproductive_ages") %in% names(demography))) {
    stop("Prepare `demography` with cohort_demography_from_spd().", call. = FALSE)
  }
  first <- min(demography$years)
  last <- max(demography$years)
  minimum <- first + n * max(demography$reproductive_ages$max_age)
  if (minimum > last) {
    stop("Insufficient simulation history for depth ", n, ": the conservative minimum sample cohort is ",
      minimum, ", after the latest census year ", last,
      ". Run a longer simulation or reduce the pedigree depth.", call. = FALSE)
  }
  for (name in c("cA", "cB")) {
    cohort <- values[[name]]
    if (cohort < minimum || cohort > last) {
      stop("`", name, "` = ", cohort, " is outside the supported sample-cohort range ", minimum,
        "--", last, " for depth ", n, ". The conservative minimum cohort is ", minimum,
        "; move the sample forward if too early, or extend the simulation if too late.", call. = FALSE)
    }
  }
  results <- nodes <- list()
  for (tree in c("A", "B")) {
    root <- if (tree == "A") cA else cB
    distributions <- vector("list", 2^(n + 1) - 1)
    for (index in seq_along(distributions)) {
      generation <- floor(log2(index))
      sex <- if (index == 1) NA_character_ else if (index %% 2 == 0) "M" else "F"
      if (index == 1) {
        dist <- data.frame(cohort = root, probability = 1)
      } else {
        child <- distributions[[floor(index / 2)]]
        ages <- demography$reproductive_ages[demography$reproductive_ages$sex == sex, ]
        cohorts <- seq.int(min(child$cohort) - ages$max_age, max(child$cohort) - ages$min_age)
        dist <- data.frame(cohort = cohorts, probability = 0)
        for (i in which(child$probability > 0)) {
          year <- child$cohort[i]
          probs <- demography$E[sex, as.character(year), ]
          if (any(!is.finite(probs)) || abs(sum(probs) - 1) > 1e-10) {
            stop("Undefined reproductive distribution for sex ", sex, " in year ", year,
              " while tracing ", tree, "-", index, ". Check reproductive abundance and schedules.", call. = FALSE)
          }
          support <- which(probs > 0)
          parent_cohorts <- year - as.numeric(dimnames(demography$E)$age[support])
          position <- match(parent_cohorts, dist$cohort)
          dist$probability[position] <- dist$probability[position] + child$probability[i] * probs[support]
        }
      }
      if (abs(sum(dist$probability) - 1) > 1e-10) stop("Cohort probability mass was lost.", call. = FALSE)
      distributions[[index]] <- dist
      individual <- if (index == 1) tree else paste0(tree, "-", index)
      node <- data.frame(tree = tree, individual = individual, index = index,
        generation = generation, sex = sex, cohort_min = min(dist$cohort), cohort_max = max(dist$cohort))
      nodes[[length(nodes) + 1L]] <- node
      results[[length(results) + 1L]] <- cbind(node[rep(1, nrow(dist)), 1:5], dist)
    }
  }
  nodes <- do.call(rbind, nodes)
  probabilities <- do.call(rbind, results)
  rownames(nodes) <- rownames(probabilities) <- NULL
  list(nodes = nodes, cohort_probabilities = probabilities, demography = demography)
}
