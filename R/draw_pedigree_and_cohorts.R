#' Draw two pedigrees and their possible ancestral cohorts
#'
#' Male ancestors are squares on the left of each child; females are circles
#' on the right. Samples are labelled A and B, with other slots numbered using
#' father = 2i and mother = 2i+1. The oldest generation's top aligns with the
#' upper edge of the earliest reference cohort. Pedigree generation heights
#' are illustrative; cohort-bar cells are aligned to actual birth years.
#'
#' @inheritParams layout_pedigree_and_cohorts
#' @inheritParams highlight_pedigree_cohorts
#' @param CohortBoxTextPad Fraction of each cohort cell reserved on each edge.
#' @param showCohorts Overlay possible cohorts on a pale reference pedigree.
#' Defaults to TRUE if a focal individual, sample_times, or shared_H11 is specified,
#' otherwise FALSE.
#' @param focalFill,reachableFill Background colours for focal cohort cells and
#' their reachable ancestral cells. Focal numbers are black for contrast.
#' @param showGrid Draw horizontal grid lines at cohort multiples of five,
#' with stronger lines at multiples of ten. Axes are always hidden.
#' @param maleColour,femaleColour Base colours used when generationPalette is NULL.
#' @param generationFade Target increase in white per generation beyond the
#' parents (0 to 1). If the oldest generation would exceed 0.5, the increments
#' are reduced uniformly to retain contrast and distinct generation colours.
#' Samples and parents use the base sex colours. Set to 0 for uniform colours.
#' Only used when generationPalette is NULL.
#' @param generationPalette Diverging ColorBrewer palette name. The default
#' "RdBu" supplies 2*n colours: red for females and blue for males, with
#' parents at the outer ends and older generations nearer the center. Samples
#' use the parental colours. Two endpoint colours are used when n is 0 or 1;
#' palettes are interpolated if 2*n exceeds the supported number of colours.
#' NULL selects the previous maleColour/femaleColour fading scheme.
#' @param cohortBarLinewidth Width of outlines and separators for all cohort
#' bars, including the black reference bar. Default 0.18 is 20 percent thicker
#' than the original 0.15.
#' @param cohortBarBold Use bold numbers in all cohort bars.
#' @param cohortBarBlackText Use black cohort numbers while retaining the
#' sex- and generation-specific bar colours. Defaults to FALSE.
#' @param sample_times Optional data frame or tibble with columns `indiv`,
#' `cohort`, and `age`. Each row marks sampling at birth cohort plus age.
#' Individuals use the same names as focal selections. Ages must be nonnegative
#' integers, and birth cohorts must be present in the individual's cohort bar.
#' Supplying this table displays cohort bars by default. Arrows run from the
#' birth-cell center to one bar-width beyond the outside edge of the bar,
#' left for A's pedigree or right for B's,
#' ending at the sampling year's vertical position.
#' @param samplingArrowLinewidth Sampling-arrow line width, defaulting to
#' 1.4 times cohortBarLinewidth.
#' @param samplingArrowColour Sampling-arrow colour, defaulting to orange.
#' @param shared_H11 Optional data frame or tibble with `focalA`, `focalB`,
#' `cohortA`, `cohortB`, `ancestorA`, and `ancestorB` columns. Ancestors must be
#' same-sex immediate parents of their respective focal individuals. The optional
#' list column `ancestor_cohorts` selects candidate birth cohorts; omit it, or
#' use NULL for a row, to select all mutually feasible cohorts. Colored solid
#' lines show alternative descent paths, matching dotted lines identify shared
#' ancestors, and parallel colored lines retain each alternative through the
#' survival interval. Black junction rectangles are one-third of a cohort cell
#' high and one-quarter of a bar wide. Lines meet these junctions at equally
#' spaced heights, from one-sixth of a cell above center to one-sixth below;
#' a single candidate uses the center.
#' This is a schematic of possibilities, not a confirmed relationship. Focal
#' shading is controlled separately. Supplying this table shows cohort bars by default.
#' @param sharedH11Linewidth Width of shared-ancestor and survival lines.
#' @param showGuides Draw thin dashed vertical guides through each pedigree
#' member, underneath all other layers.
#' @param pedigreeColour Colour of pedigree outlines, labels and edges. NULL
#' uses black, or grey70 when cohort bars are overlaid.
#' @return A ggplot object. Use [layout_pedigree_and_cohorts()] to obtain the
#' underlying coordinates and feasible cohort intervals separately.
#' @examples
#' draw_pedigree_and_cohorts()
#' draw_pedigree_and_cohorts(showCohorts = TRUE)
#' draw_pedigree_and_cohorts(focal_indiv_A = "A-2", focal_cohort_A = 51,
#'   focal_indiv_B = "B-2", focal_cohort_B = 54)
#' shared <- data.frame(focalA = "A", focalB = "B", cohortA = 56, cohortB = 60,
#'   ancestorA = "A-2", ancestorB = "B-2")
#' draw_pedigree_and_cohorts(shared_H11 = shared, cohortBarBlackText = TRUE)
#' @md
#' @export
draw_pedigree_and_cohorts <- function(
  Am = 10, Af = 12, Mm = 3, Mf = 6,
  cA = 56, cB = 60, aA = 5, aB = 3, n = 3,
  xPosBar = 0, CohortBoxWidth = 1, CohortBoxTextPad = 0.05,
  xPosA = NULL, xPosB = NULL, xPedSep = 16,
  sexA = "F", sexB = "F", nodeWidth = NULL, nodeHeight = NULL,
  showCohorts = !is.null(focal_indiv_A) || !is.null(focal_indiv_B) || !is.null(sample_times) || !is.null(shared_H11), showGrid = TRUE,
  maleColour = "#2166AC", femaleColour = "#B2182B", pedigreeColour = NULL,
  SampleSepFactor = 2, generationFade = 0.18, showGuides = TRUE,
  cohortBarLinewidth = 0.18, cohortBarBold = FALSE, generationPalette = "RdBu",
  focal_indiv_A = NULL, focal_cohort_A = NULL,
  focal_indiv_B = NULL, focal_cohort_B = NULL,
  focalFill = "#FF00FF", reachableFill = "grey88",
  cohortBarBlackText = FALSE, sample_times = NULL,
  samplingArrowLinewidth = 1.4 * cohortBarLinewidth, samplingArrowColour = "orange",
  shared_H11 = NULL, sharedH11Linewidth = 1.4 * cohortBarLinewidth
) {
  if (!is.numeric(generationFade) || length(generationFade) != 1L ||
      !is.finite(generationFade) || generationFade < 0 || generationFade > 1) {
    stop("`generationFade` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(cohortBarBold) || length(cohortBarBold) != 1L || is.na(cohortBarBold)) {
    stop("`cohortBarBold` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(cohortBarBlackText) || length(cohortBarBlackText) != 1L ||
      is.na(cohortBarBlackText)) {
    stop("`cohortBarBlackText` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(samplingArrowLinewidth) || length(samplingArrowLinewidth) != 1L ||
      !is.finite(samplingArrowLinewidth) || samplingArrowLinewidth <= 0) {
    stop("`samplingArrowLinewidth` must be one positive finite number.", call. = FALSE)
  }
  palette_colours <- NULL
  if (!is.null(generationPalette)) {
    info <- RColorBrewer::brewer.pal.info
    if (length(generationPalette) != 1L || !is.character(generationPalette) ||
        is.na(generationPalette) || !generationPalette %in% rownames(info) ||
        info[generationPalette, "category"] != "div") {
      stop("`generationPalette` must be a diverging ColorBrewer palette name or NULL.", call. = FALSE)
    }
  }
  layout <- layout_pedigree_and_cohorts(Am = Am, Af = Af, Mm = Mm, Mf = Mf,
    cA = cA, cB = cB, aA = aA, aB = aB, n = n, xPosBar = xPosBar,
    xPosA = xPosA, xPosB = xPosB, xPedSep = xPedSep,
    CohortBoxWidth = CohortBoxWidth, sexA = sexA, sexB = sexB,
    nodeWidth = nodeWidth, nodeHeight = nodeHeight, SampleSepFactor = SampleSepFactor)
  layout <- highlight_pedigree_cohorts(layout,
    focal_indiv_A = focal_indiv_A, focal_cohort_A = focal_cohort_A,
    focal_indiv_B = focal_indiv_B, focal_cohort_B = focal_cohort_B)
  sampling_arrows <- pedigree_sampling_arrows(layout, sample_times)
  if (!is.numeric(sharedH11Linewidth) || length(sharedH11Linewidth) != 1L ||
      !is.finite(sharedH11Linewidth) || sharedH11Linewidth <= 0) {
    stop("`sharedH11Linewidth` must be one positive finite number.", call. = FALSE)
  }
  shared_relationships <- pedigree_shared_H11(layout, shared_H11)
  for (relation in shared_relationships) {
    layout$xlim <- range(c(layout$xlim, relation$points$x - CohortBoxWidth / 4,
      relation$points$x + CohortBoxWidth / 4))
  }
  if (nrow(sampling_arrows)) {
    latest <- max(sampling_arrows$sampling_cohort)
    if (latest > layout$BarLate) {
      extra <- layout$reference_bar[rep(1L, latest - layout$BarLate), ]
      extra$cohort <- seq.int(layout$BarLate + 1, latest)
      extra$ymid <- layout$BarLate - extra$cohort
      layout$reference_bar <- rbind(layout$reference_bar, extra)
      grid <- data.frame(cohort = extra$cohort[extra$cohort %% 5 == 0])
      grid$y <- layout$BarLate - grid$cohort
      grid$major <- grid$cohort %% 10 == 0
      layout$grid <- rbind(layout$grid, grid)
      layout$ylim[1] <- min(layout$ylim[1], layout$BarLate - latest - 1)
    }
    layout$xlim <- range(c(layout$xlim, sampling_arrows$xend))
  }
  if (is.null(pedigreeColour)) {
    pedigreeColour <- if (showCohorts) "grey70" else "black"
  }
  if (!is.null(generationPalette)) {
    size <- 2 * max(1, n)
    maximum <- info[generationPalette, "maxcolors"]
    palette_colours <- RColorBrewer::brewer.pal(max(3, min(size, maximum)), generationPalette)
    if (size == 2) palette_colours <- palette_colours[c(1, length(palette_colours))]
    if (size > maximum) palette_colours <- grDevices::colorRampPalette(palette_colours)(size)
  }
  fontface <- if (cohortBarBold) "bold" else "plain"
  plot <- ggplot2::ggplot()
  if (showGuides) plot <- plot + pedigree_guide_layer(layout, pedigreeColour)
  if (showGrid) {
    plot <- plot + list(
      ggplot2::geom_hline(data = layout$grid[!layout$grid$major, , drop = FALSE],
        ggplot2::aes(yintercept = y), colour = "grey94", linewidth = 0.2),
      ggplot2::geom_hline(data = layout$grid[layout$grid$major, , drop = FALSE],
        ggplot2::aes(yintercept = y), colour = "grey86", linewidth = 0.3)
    )
  }
  plot <- plot + draw_layer_pedigree(layout, colour = pedigreeColour, showGuides = FALSE)
  if (showCohorts) {
    for (sex in c("M", "F")) {
      base_colour <- if (sex == "M") maleColour else femaleColour
      for (generation in sort(unique(layout$cohort_bars$generation))) {
        bars <- layout$cohort_bars[layout$cohort_bars$sex == sex &
          layout$cohort_bars$generation == generation, ]
        if (is.null(palette_colours)) {
          oldest_fade <- min(0.5, max(0, n - 1) * generationFade)
          fade <- oldest_fade * max(0, generation - 1) / max(1, n - 1)
          rgb <- grDevices::colorRamp(c(base_colour, "white"))(fade)
          colour <- grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
        } else {
          index <- max(1, generation)
          if (sex == "M") index <- length(palette_colours) - index + 1
          colour <- palette_colours[index]
        }
        bars$cell_fill <- ifelse(bars$focal, focalFill,
          ifelse(bars$reachable, reachableFill, NA_character_))
        bars$cell_text_colour <- ifelse(bars$focal | cohortBarBlackText, "black", colour)
        plot <- plot + draw_layer_cohort_bar(bars, pad = CohortBoxTextPad,
          colour = colour, fill = "white", group = "id",
          linewidth = cohortBarLinewidth, fontface = fontface,
          cellFill = "cell_fill", textColour = "cell_text_colour")
      }
    }
  }
  plot <- plot + draw_layer_cohort_bar(layout$reference_bar, pad = CohortBoxTextPad,
    fill = "white", linewidth = cohortBarLinewidth, fontface = fontface)
  if (nrow(sampling_arrows)) {
    plot <- plot + ggplot2::geom_segment(data = sampling_arrows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = samplingArrowColour, linewidth = samplingArrowLinewidth,
      arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed"),
      inherit.aes = FALSE)
  }
  plot <- plot + pedigree_shared_H11_layers(layout, shared_relationships, sharedH11Linewidth)
  plot +
    ggplot2::coord_fixed(ratio = layout$aspect_ratio, xlim = layout$xlim,
      ylim = layout$ylim, expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 8, 8, 8))
}
