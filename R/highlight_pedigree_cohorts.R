#' Mark focal cohorts and reachable ancestral cohorts
#'
#' Conditions an ancestral slot on one exact birth cohort, then propagates
#' possible birth years to its ancestors using sex-specific reproductive ages.
#' Other branches and the other pedigree are unaffected. Reachability describes
#' possibilities, not probabilities or identities of shared individuals.
#'
#' @param layout A list returned by [layout_pedigree_and_cohorts()].
#' @param focal_indiv_A,focal_indiv_B Focal individual in each pedigree. Use
#' "A" or "B" for a sample, or a numbered slot such as "A-2" (father) or
#' "B-6" (maternal grandfather). Explicit sample indices "A-1" and "B-1"
#' are also accepted. NULL leaves that pedigree unselected.
#' @param focal_cohort_A,focal_cohort_B Exact integer birth cohort of the focal
#' individual. Must be an existing cell in that individual's cohort bar.
#' Supply an individual and cohort together, or leave both NULL.
#' @return The layout with `focal`, `reachable`, and `highlight` columns added
#' to `cohort_bars`, and a `focal_selections` table. `reachable` includes the
#' focal cell; `highlight` is "focal", "reachable", or "none". Previous
#' selections are replaced. All coordinates and full cohort ranges are retained.
#' @examples
#' layout <- layout_pedigree_and_cohorts()
#' layout <- highlight_pedigree_cohorts(layout,
#'   focal_indiv_A = "A-2", focal_cohort_A = 51)
#' subset(layout$cohort_bars, reachable)
#' @md
#' @export
highlight_pedigree_cohorts <- function(layout,
  focal_indiv_A = NULL, focal_cohort_A = NULL,
  focal_indiv_B = NULL, focal_cohort_B = NULL
) {
  bars <- layout$cohort_bars
  nodes <- layout$nodes
  bars$focal <- FALSE
  bars$reachable <- FALSE
  selections <- data.frame(tree = character(), individual = character(),
    id = character(), cohort = numeric(), stringsAsFactors = FALSE)
  for (tree in c("A", "B")) {
    individual <- if (tree == "A") focal_indiv_A else focal_indiv_B
    cohort <- if (tree == "A") focal_cohort_A else focal_cohort_B
    if (is.null(individual) && is.null(cohort)) next
    if (is.null(individual) || is.null(cohort)) {
      stop("Supply both focal_indiv_", tree, " and focal_cohort_", tree, ".", call. = FALSE)
    }
    if (!is.character(individual) || length(individual) != 1L || is.na(individual) ||
        !grepl(paste0("^", tree, "(-[1-9][0-9]*)?$"), individual)) {
      stop("Use '", tree, "' or '", tree, "-<index>' for focal_indiv_", tree, ".", call. = FALSE)
    }
    if (!is.numeric(cohort) || length(cohort) != 1L || !is.finite(cohort) || cohort != floor(cohort)) {
      stop("focal_cohort_", tree, " must be one finite integer.", call. = FALSE)
    }
    index <- if (individual == tree) 1 else as.numeric(sub(paste0("^", tree, "-"), "", individual))
    selected <- which(nodes$tree == tree & nodes$index == index)
    if (!length(selected)) {
      stop("'", individual, "' is not present at the requested pedigree depth.", call. = FALSE)
    }
    id <- nodes$id[selected]
    cell <- bars$id == id & bars$cohort == cohort
    if (!any(cell)) {
      stop("Cohort ", cohort, " is not possible for ", individual,
        "; choose a cohort from ", nodes$cohort_min[selected], " to ",
        nodes$cohort_max[selected], ".", call. = FALSE)
    }
    low <- high <- rep(NA_real_, nrow(nodes))
    low[selected] <- high[selected] <- cohort
    # The layout's child index precedes both parental indices.
    for (i in which(nodes$tree == tree & nodes$index > index)) {
      child <- match(paste(tree, nodes$parent_index[i], sep = ":"), nodes$id)
      if (is.na(child) || is.na(low[child])) next
      age <- layout$reproductive_ages[layout$reproductive_ages$sex == nodes$sex[i], ]
      low[i] <- max(nodes$cohort_min[i], low[child] - age$max_age)
      high[i] <- min(nodes$cohort_max[i], high[child] - age$min_age)
    }
    node <- match(bars$id, nodes$id)
    reachable <- !is.na(low[node]) & bars$cohort >= low[node] & bars$cohort <= high[node]
    bars$reachable <- bars$reachable | reachable
    bars$focal <- bars$focal | cell
    selections <- rbind(selections, data.frame(tree = tree, individual = individual,
      id = id, cohort = cohort, stringsAsFactors = FALSE))
  }
  bars$highlight <- ifelse(bars$focal, "focal", ifelse(bars$reachable, "reachable", "none"))
  layout$cohort_bars <- bars
  layout$focal_selections <- selections
  layout
}
