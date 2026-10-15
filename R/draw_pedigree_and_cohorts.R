


# helper function to plot a cohort bar
draw_layer_cohort_bar <- function(BT, pad) {
  BT_rect <- BT %>%
    summarise(
      xmin = min(xlo),
      xmax =  max(xhi),
      ymin = min(ymid - 0.5),
      ymax = max(ymid + 0.5)
    )
  list(
    geom_rect(
    data = BT_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    colour = "black",
    linewidth = 0.1
    ),
    geom_segment(
      data = BT,
      aes(x = xlo, xend = xhi, y = ymid - 0.5),
      colour = "black",
      linewidth = 0.1
    ),
    ggfittext::geom_fit_text(
      data = BT,
      aes(
        xmin = xlo + pad * (xhi - xlo),
        xmax = xhi - pad * (xhi - xlo),
        ymin = ymid - 0.5 + pad,
        ymax = ymid + 0.5 - pad,
        label = cohort
      ),
      min.size = 0,
      grow = TRUE,
      padding.x = grid::unit(0, "mm"),
      padding.y = grid::unit(0, "mm")
    )
  )
}




#' Draw pedigrees and the associated cohorts different members could belong to
#'
#' This provides a retrospective view of which cohorts the ancestors of two
#' distinct individuals (belonging to particular cohorts) could belong to. The
#' hope is that this will help provide some insight into how we might efficiently
#' design some algorithms to properly calculate multi-generational kinship
#' probabilities in age structured populations.  Note that we have to provide it
#' with some information about the life history of the population.  Most notably,
#' we are going to need max reproductive ages and minimum ages at maturity.  Because we are going
#' to be keeping track of the sex of ancestors, and in order to be really general
#' about it, we will assume that both max and min reproductive ages are sex specific.
#' The defaults are set for a population where males mature by 3, females by 6, but
#' males are dead/non-reproductive by 10 and females are done/non-reproductive by 12.
#' @param Am maximum reproductive age of males
#' @param Af maximum reproductive age of females
#' @param Mm minimum reproductive age of males
#' @param Mf minimum reproductive age of females
#' @param cA cohort of sample A
#' @param cB cohort of sample B
#' @param aA age at which A is sampled
#' @param aB age at which B is sampled
#' @param n number of generations back to consider
#' @param xPosBar center x-position of the cohort "scale bar"
#' @param xPosA x-position for sample A
#' @param xPosB x-position for sample B
#' @param xPedSep the amount (positive for mother and negative for father)
#' that the father and mother are moved in their x-dimension relative to the
#' sample.  Note that for subsequent generations the change in x is 1/2 of this
#' for every generation further back we go.
#' @param CohortBoxWidth total width of the cohort reference bar and the pedigree
#' ancestor cohort bars
#' @param CohortBoxTextPad fraction of the cohort bar box to leave as margins for the
#' text/numerals inside
#' @export
draw_pedigree_and_cohorts <- function(
  Am = 10,
  Af = 12,
  Mm = 3,
  Mf = 6,
  cA = 56,
  cB = 60,
  aA = 5,
  aB = 3,
  n = 3,
  xPosBar = 0,
  CohortBoxWidth = 1,
  CohortBoxTextPad = 0.05
) {

  # first, let us calculate the y-limits on the cohort bar that we will want.
  # These are for reference and will go from the earliest possible ancestor to
  # the latest time of sampling.
  BarEarly <- min(cA, cB) - n * max(Am, Af)
  BarLate <- max(cA + aA, cB + aB)

  # now, we need a function for translating cohort to height on the y-axis. Time
  # runs from the top down, so things are reversed. We will let BarLate correspond
  # to a y of 0, and everything before that is positive y.
  # ZeroCohort here is the cohort we want at y=0.  C is a vector of cohort times)
  getCohortY <- function(ZeroCohort = BarLate, C) {
    -(C - ZeroCohort)
  }


  # Make a tibble with needed parameters for drawing the cohort reference bar
  CohortBarTib <- tibble(cohort = BarEarly:BarLate) %>%
    mutate(
      ymid = getCohortY(C = cohort),
      xmid = xPosBar,
      xwidth = CohortBoxWidth,
      xlo = xmid - xwidth / 2,
      xhi = xmid + xwidth / 2,
    )

  ggplot() +
    draw_layer_cohort_bar(CohortBarTib, pad = CohortBoxTextPad) +
    xlim(-1,40)

}

