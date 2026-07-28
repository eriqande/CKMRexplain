# Lightweight CKMR model graph -------------------------------------------------
#
# The functions in this file define the first "grammar object" for CKMRexplain.
# The object is deliberately plain: two data frames, one for nodes and one for
# edges, plus a small list for richer node-level data such as distributions.
# That makes the model easy to print, inspect, test, and eventually hand to
# plotting, simulation, Shiny, or TMB helper functions.

ckmr_model <- function(nodes = NULL, edges = NULL, node_data = list()) {
  if (is.null(nodes)) {
    nodes <- empty_node_table()
  }
  if (is.null(edges)) {
    edges <- empty_edge_table()
  }

  structure(
    list(
      nodes = nodes,
      edges = edges,
      node_data = node_data
    ),
    class = "ckmr_model"
  )
}

empty_node_table <- function() {
  data.frame(
    node_id = character(),
    node_type = character(),
    label = character(),
    sex = character(),
    stage = character(),
    age = numeric(),
    time = numeric(),
    abundance = numeric(),
    fecundity = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_edge_table <- function() {
  data.frame(
    edge_id = character(),
    from = character(),
    to = character(),
    edge_type = character(),
    label = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )
}

print.ckmr_model <- function(x, ...) {
  cat("<ckmr_model>\n")
  cat("  nodes:", nrow(x$nodes), "\n")
  cat("  edges:", nrow(x$edges), "\n")

  if (nrow(x$nodes) > 0) {
    cat("\nNode types:\n")
    print(table(x$nodes$node_type))
  }

  invisible(x)
}

add_population_node <- function(model,
                                sex = NA_character_,
                                stage,
                                time = 0,
                                ages = NA_real_,
                                abundance = NA_real_,
                                fecundity = NA_real_,
                                label = NULL,
                                id = NULL,
                                covariates = list()) {
  stopifnot_ckmr_model(model)

  # We expand over every requested time-age combination. This is intentionally
  # simple, but it means Leslie-matrix-like outputs can already be dropped in as
  # matrices with rows = times and columns = ages.
  grid <- expand.grid(
    time = time,
    age = ages,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  abundance <- expand_value(abundance, nrow(grid), "abundance")
  fecundity <- expand_value(fecundity, nrow(grid), "fecundity")

  if (is.null(id)) {
    node_id <- make_population_ids(sex, stage, grid$time, grid$age)
  } else if (length(id) == 1 && nrow(grid) > 1) {
    node_id <- paste(id, seq_len(nrow(grid)), sep = "_")
  } else {
    node_id <- expand_value(id, nrow(grid), "id")
  }

  if (is.null(label)) {
    label <- make_population_labels(sex, stage, grid$time, grid$age)
  } else {
    label <- expand_value(label, nrow(grid), "label")
  }

  new_nodes <- data.frame(
    node_id = node_id,
    node_type = "population",
    label = label,
    sex = rep(as.character(sex), nrow(grid)),
    stage = rep(as.character(stage), nrow(grid)),
    age = grid$age,
    time = grid$time,
    abundance = abundance,
    fecundity = fecundity,
    stringsAsFactors = FALSE
  )

  model$nodes <- rbind(model$nodes, new_nodes)

  # Rich node data are kept outside the flat node table. For now this is mostly
  # for covariate distributions, such as adult length distributions that affect
  # fecundity or sampling selectivity.
  for (node in node_id) {
    model$node_data[[node]] <- list(covariates = covariates)
  }

  model
}

add_process_node <- function(model,
                             id,
                             node_type,
                             label = id,
                             time = NA_real_) {
  stopifnot_ckmr_model(model)

  if (id %in% model$nodes$node_id) {
    return(model)
  }

  new_node <- data.frame(
    node_id = id,
    node_type = node_type,
    label = label,
    sex = NA_character_,
    stage = NA_character_,
    age = NA_real_,
    time = time,
    abundance = NA_real_,
    fecundity = NA_real_,
    stringsAsFactors = FALSE
  )

  model$nodes <- rbind(model$nodes, new_node)
  model$node_data[[id]] <- list(covariates = list())
  model
}

add_reproduction <- function(model,
                             from,
                             to,
                             reproduction_node = "mating",
                             label = "mating",
                             weight = NA_real_) {
  stopifnot_ckmr_model(model)

  model <- add_process_node(
    model,
    id = reproduction_node,
    node_type = "reproduction",
    label = label
  )

  input_edges <- data.frame(
    edge_id = paste("repro_in", from, reproduction_node, sep = "_"),
    from = from,
    to = reproduction_node,
    edge_type = "reproduction_input",
    label = "",
    weight = weight,
    stringsAsFactors = FALSE
  )

  output_edge <- data.frame(
    edge_id = paste("birth_out", reproduction_node, to, sep = "_"),
    from = reproduction_node,
    to = to,
    edge_type = "birth_output",
    label = "",
    weight = weight,
    stringsAsFactors = FALSE
  )

  model$edges <- rbind(model$edges, input_edges, output_edge)
  model
}

add_sampling <- function(model,
                         from,
                         sample = NULL,
                         label = "sample",
                         fraction = NA_real_) {
  stopifnot_ckmr_model(model)

  if (is.null(sample)) {
    sample <- paste0(from, "_sample")
  }

  model <- add_process_node(
    model,
    id = sample,
    node_type = "sample",
    label = label
  )

  new_edge <- data.frame(
    edge_id = paste("sample", from, sample, sep = "_"),
    from = from,
    to = sample,
    edge_type = "sampling",
    label = "",
    weight = fraction,
    stringsAsFactors = FALSE
  )

  model$edges <- rbind(model$edges, new_edge)
  model
}

add_survival <- function(model,
                         from,
                         to,
                         label = "",
                         survival = NA_real_) {
  stopifnot_ckmr_model(model)

  new_edges <- data.frame(
    edge_id = paste("survival", from, to, sep = "_"),
    from = from,
    to = to,
    edge_type = "survival",
    label = label,
    weight = survival,
    stringsAsFactors = FALSE
  )

  model$edges <- rbind(model$edges, new_edges)
  model
}

add_recruitment <- function(model,
                            from,
                            to,
                            label = "",
                            recruitment = NA_real_) {
  stopifnot_ckmr_model(model)

  new_edges <- data.frame(
    edge_id = paste("recruitment", from, to, sep = "_"),
    from = from,
    to = to,
    edge_type = "recruitment",
    label = label,
    weight = recruitment,
    stringsAsFactors = FALSE
  )

  model$edges <- rbind(model$edges, new_edges)
  model
}

demo_two_stage_po <- function(include_size = TRUE) {
  female_lengths <- NULL
  male_lengths <- NULL

  if (include_size) {
    female_lengths <- list(length = normal_covariate(mean = 82, sd = 9))
    male_lengths <- list(length = normal_covariate(mean = 88, sd = 11))
  }

  ckmr_model() |>
    add_population_node(
      id = "adult_F",
      sex = "F",
      stage = "adult",
      abundance = 1000,
      fecundity = 1.35,
      label = "adult females",
      covariates = female_lengths
    ) |>
    add_population_node(
      id = "adult_M",
      sex = "M",
      stage = "adult",
      abundance = 800,
      fecundity = 1.00,
      label = "adult males",
      covariates = male_lengths
    ) |>
    add_process_node(
      id = "offspring_pool",
      node_type = "offspring_pool",
      label = "offspring"
    ) |>
    add_reproduction(
      from = c("adult_F", "adult_M"),
      to = "offspring_pool",
      reproduction_node = "mating",
      label = "mating"
    ) |>
    add_sampling(from = "adult_F", sample = "adult_F_sample") |>
    add_sampling(from = "adult_M", sample = "adult_M_sample") |>
    add_sampling(from = "offspring_pool", sample = "offspring_sample")
}

demo_age_time_model <- function(times = 0:2,
                                ages = 1:8,
                                include_samples = TRUE) {
  model <- ckmr_model()

  # Give abundance a gentle deterministic pattern so the abundance overlay has
  # something to show. Rows are times and columns are ages, matching the shape
  # that Leslie-matrix simulation output would often have.
  abundance_surface <- outer(
    times,
    ages,
    function(tt, aa) 1000 * exp(-0.16 * aa) * (1 + 0.04 * tt)
  )
  fecundity_surface <- outer(
    times,
    ages,
    function(tt, aa) aa / max(ages)
  )

  for (sex in c("F", "M")) {
    model <- add_population_node(
      model,
      sex = sex,
      stage = "adult",
      time = times,
      ages = ages,
      abundance = abundance_surface,
      fecundity = if (sex == "F") fecundity_surface else 1,
      label = NULL
    )
  }

  for (tt in times) {
    repro_node <- paste0("mating_t", tt)
    offspring_node <- paste0("offspring_t", tt)

    model <- add_process_node(
      model,
      id = repro_node,
      node_type = "reproduction",
      label = "mate",
      time = tt
    )
    model <- add_process_node(
      model,
      id = offspring_node,
      node_type = "offspring_pool",
      label = "offs.",
      time = tt
    )

    female_nodes <- make_population_ids("F", "adult", tt, ages)
    male_nodes <- make_population_ids("M", "adult", tt, ages)

    model <- add_reproduction(
      model,
      from = c(female_nodes, male_nodes),
      to = offspring_node,
      reproduction_node = repro_node,
      label = "mate"
    )
  }

  if (length(times) > 1) {
    for (tt in times[-length(times)]) {
      offspring_node <- paste0("offspring_t", tt)

      for (sex in c("F", "M")) {
        model <- add_recruitment(
          model,
          from = offspring_node,
          to = make_population_ids(sex, "adult", tt + 1, min(ages)),
          label = "R"
        )
      }

      for (sex in c("F", "M")) {
        for (age in ages[-length(ages)]) {
          model <- add_survival(
            model,
            from = make_population_ids(sex, "adult", tt, age),
            to = make_population_ids(sex, "adult", tt + 1, age + 1),
            label = "S"
          )
        }
      }
    }
  }

  if (include_samples) {
    last_time <- max(times)
    sample_ages <- c(min(ages), max(ages))
    for (sex in c("F", "M")) {
      for (age in sample_ages) {
        node <- make_population_ids(sex, "adult", last_time, age)
        model <- add_sampling(
          model,
          from = node,
          sample = paste0(node, "_sample")
        )
      }
    }
    model <- add_sampling(
      model,
      from = paste0("offspring_t", last_time),
      sample = paste0("offspring_t", last_time, "_sample")
    )
  }

  model
}

normal_covariate <- function(mean, sd, n = 80) {
  # Store a deterministic curve, not random draws, so examples and tests are
  # reproducible. The plotting layer can draw this as a tiny density curve.
  x <- seq(mean - 3 * sd, mean + 3 * sd, length.out = n)
  density <- stats::dnorm(x, mean = mean, sd = sd)
  data.frame(x = x, density = density, stringsAsFactors = FALSE)
}

stopifnot_ckmr_model <- function(model) {
  if (!inherits(model, "ckmr_model")) {
    stop("`model` must be a <ckmr_model> object.", call. = FALSE)
  }
}

expand_value <- function(value, n, arg) {
  if (length(value) == 1) {
    return(rep(value, n))
  }

  if (is.matrix(value) || is.data.frame(value)) {
    value <- as.vector(as.matrix(value))
  }

  if (length(value) != n) {
    stop(
      "`", arg, "` must have length 1 or length ", n,
      " after expansion.",
      call. = FALSE
    )
  }

  value
}

make_population_ids <- function(sex, stage, time, age) {
  sex_part <- ifelse(is.na(sex), "all", sex)
  age_part <- ifelse(is.na(age), "na", paste0("a", age))
  paste(sex_part, stage, paste0("t", time), age_part, sep = "_")
}

make_population_labels <- function(sex, stage, time, age) {
  sex_part <- ifelse(is.na(sex), "", paste0(sex, " "))
  age_part <- ifelse(is.na(age), "", paste0(", age ", age))
  paste0(sex_part, stage, "\nt=", time, age_part)
}
