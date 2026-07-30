\name{ckmr_model}
\alias{ckmr_model}
\alias{add_population_node}
\alias{add_process_node}
\alias{add_recruitment}
\alias{add_reproduction}
\alias{add_sampling}
\alias{add_survival}
\alias{demo_age_time_model}
\alias{demo_half_thiatic_path_model}
\alias{demo_half_thiatic_younger_uncle_path_model}
\alias{demo_mhsp_path_model}
\alias{demo_two_stage_po}
\alias{layout_ckmr_model}
\alias{normal_covariate}
\alias{plot_ckmr_age_time}
\alias{plot_ckmr_model}
\title{Build and plot lightweight CKMR model graphs}
\usage{
ckmr_model(nodes = NULL, edges = NULL, node_data = list())

add_population_node(
  model,
  sex = NA_character_,
  stage,
  time = 0,
  ages = NA_real_,
  abundance = NA_real_,
  fecundity = NA_real_,
  label = NULL,
  id = NULL,
  covariates = list()
)

add_process_node(model, id, node_type, label = id, time = NA_real_)

add_recruitment(model, from, to, label = "", recruitment = NA_real_)

add_reproduction(
  model,
  from,
  to,
  reproduction_node = "mating",
  label = "mating",
  weight = NA_real_
)

add_sampling(model, from, sample = NULL, label = "sample", fraction = NA_real_)

add_survival(model, from, to, label = "", survival = NA_real_)

demo_age_time_model(times = 0:2, ages = 1:8, include_samples = TRUE)

demo_half_thiatic_path_model(parent_birth_times = 12:14)

demo_half_thiatic_younger_uncle_path_model(parent_birth_times = 9:11)

demo_mhsp_path_model()

demo_two_stage_po(include_size = TRUE)

layout_ckmr_model(model, layout = "two_stage_po")

normal_covariate(mean, sd, n = 80)

plot_ckmr_age_time(
  model,
  overlay = c("none", "abundance", "fecundity"),
  show_reproduction_inputs = TRUE,
  highlight_edges = attr(model, "highlight_edges", exact = TRUE),
  path_highlights = attr(model, "path_highlights", exact = TRUE),
  selected_paths = attr(model, "selected_paths", exact = TRUE)
)

plot_ckmr_model(
  model,
  overlay = c("none", "abundance", "fecundity", "length_distribution"),
  layout = "two_stage_po"
)
}
\arguments{
\item{model}{A \code{ckmr_model} object.}
\item{nodes, edges}{Optional node and edge data frames.}
\item{node_data}{A named list of richer node attributes, such as covariate distributions.}
\item{sex}{Sex label for a population node, such as \code{"F"} or \code{"M"}.}
\item{stage}{Stage label for a population node, such as \code{"adult"}.}
\item{time}{One or more time values.}
\item{ages}{One or more age values.}
\item{abundance, fecundity}{Scalar, vector, or matrix values attached to expanded population nodes.}
\item{label}{Human-readable label.}
\item{id}{Node identifier.}
\item{covariates}{Named list of covariate distributions or other node-level attributes.}
\item{node_type}{Type of process node, such as \code{"reproduction"}, \code{"sample"}, or \code{"offspring_pool"}.}
\item{from, to}{Node identifiers used to create graph edges.}
\item{sample}{Identifier for the sample node created by \code{add_sampling()}.}
\item{reproduction_node}{Identifier for the mating or reproductive process node.}
\item{weight, fraction}{Optional numeric edge weights.}
\item{recruitment}{Optional recruitment edge weight.}
\item{survival}{Optional survival edge weight.}
\item{times}{Time values for the age-time demo model.}
\item{parent_birth_times}{Intervening parent birth years to show in the half-thiatic demo.}
\item{include_samples}{Whether to include a small number of sample nodes in the age-time demo.}
\item{include_size}{Whether the demo model includes adult length distributions.}
\item{layout}{Layout name. Currently only \code{"two_stage_po"} is implemented.}
\item{mean, sd, n}{Parameters used to create a deterministic normal covariate curve.}
\item{overlay}{Visual overlay to draw on population nodes.}
\item{show_reproduction_inputs}{Whether to show all age-specific inputs to the mating boxes.}
\item{highlight_edges}{Optional edge identifiers to draw more heavily.}
\item{path_highlights}{Optional path catalog with path IDs, edge IDs, edge roles, and colors.}
\item{selected_paths}{Optional path IDs to draw from the path catalog.}
}
\description{
These functions define the first lightweight graph grammar for CKMRexplain.
The model stores nodes, edges, and richer node data separately so that the same
object can later support plotting, simulation, Shiny apps, and inference
helpers.
}
\examples{
model <- demo_two_stage_po()
model

plot_ckmr_model(model, overlay = "length_distribution")
}
