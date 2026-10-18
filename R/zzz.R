# Tell R CMD check that these names are data-frame columns used inside ggplot2
# aesthetics. This keeps the plotting code readable while avoiding noisy
# "no visible binding" notes.
utils::globalVariables(c(
  "color", "display_label", "fill_value", "label", "node_id", "sex",
  "edge_alpha", "edge_width", "linetype", "sex_label",
  "x", "x_end", "x_scaled", "x_start", "xmax", "xmin",
  "y", "y_end", "y_scaled", "y_start", "ymax", "ymin",
  "cohort", "id", "xend", "xhi", "xlo", "yend", "ymid",
  "text_xmin", "text_xmax", "text_ymin", "text_ymax"
))
