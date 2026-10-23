#' Prepare reproductive contributions from a CKMRpop simulation
#'
#' Uses the census after mortality and before reproduction. No extra survival
#' factor is applied: survivors are already represented in the census.
#' CKMRpop itself is not required to process these input objects.
#' @param spd A CKMRpop/SPIP parameter list, including max-age, sex-specific
#' prob-repro, asrf/asrp, and surv-probs entries.
#' @param census A census_postkill data frame with year, pop, age, male,
#' and female columns. Every year must explicitly record ages 1 through max-age,
#' including zero counts. Additional age-zero records are ignored.
#' @param pop Population identifier to select.
#' @return Demographic inputs with sex-by-year-by-age N and E arrays,
#' reproductive age limits, and survival schedules retained for future use.
#' E is NA for a sex/year with no reproductive contribution; attempting to use
#' that year in the pedigree calculation produces an error.
#' @export
cohort_demography_from_spd <- function(spd, census, pop = 0L) {
  if (!is.list(spd)) stop("`spd` must be a parameter list.", call. = FALSE)
  names(spd) <- gsub("\\p{Pd}", "-", names(spd), perl = TRUE)
  required <- c("max-age", "male-prob-repro", "fem-prob-repro",
    "male-asrp", "fem-asrf", "male-surv-probs", "fem-surv-probs")
  missing <- setdiff(required, names(spd))
  if (length(missing)) stop("Missing SPD entries: ", paste(missing, collapse = ", "), call. = FALSE)
  Amax <- spd[["max-age"]]
  if (!is.numeric(Amax) || length(Amax) != 1L || !is.finite(Amax) ||
      Amax < 1 || Amax != floor(Amax)) stop("`max-age` must be a positive integer.", call. = FALSE)
  for (name in required[-1]) {
    values <- spd[[name]]
    if (!is.numeric(values) || length(values) != Amax || any(!is.finite(values)) ||
        any(values < 0) || (grepl("prob|surv", name) && any(values > 1))) {
      stop("Invalid SPD schedule: ", name, ". Supply max-age finite nonnegative values; probabilities must be at most one.", call. = FALSE)
    }
  }
  columns <- c("year", "pop", "age", "male", "female")
  if (!is.data.frame(census) || !all(columns %in% names(census))) {
    stop("`census` must contain year, pop, age, male, and female columns.", call. = FALSE)
  }
  if (!is.numeric(pop) || length(pop) != 1L || !is.finite(pop) || pop != floor(pop)) {
    stop("`pop` must be one finite integer.", call. = FALSE)
  }
  for (column in columns) {
    values <- census[[column]]
    if (!is.numeric(values) || any(!is.finite(values)) ||
        (column %in% c("year", "pop", "age") && any(values != floor(values))) ||
        (column %in% c("age", "male", "female") && any(values < 0))) {
      stop("Invalid census column: ", column, ".", call. = FALSE)
    }
  }
  census <- census[census$pop == pop, columns, drop = FALSE]
  if (!nrow(census)) stop("No census records for population ", pop, ".", call. = FALSE)
  if (any(census$age > Amax)) stop("Census ages exceed SPD max-age.", call. = FALSE)
  years <- seq.int(min(census$year), max(census$year))
  census <- census[census$age >= 1, , drop = FALSE]
  keys <- paste(census$year, census$age)
  if (anyDuplicated(keys)) stop("Duplicate census year/age records for this population.", call. = FALSE)
  expected <- expand.grid(year = years, age = seq_len(Amax))
  index <- match(paste(expected$year, expected$age), keys)
  if (anyNA(index)) {
    first <- expected[which(is.na(index))[1], ]
    stop("Incomplete census: missing age ", first$age, " in year ", first$year,
      ". Supply every year and age, explicitly recording zero counts.", call. = FALSE)
  }
  N <- array(0, c(2, length(years), Amax),
    dimnames = list(sex = c("M", "F"), year = as.character(years), age = as.character(seq_len(Amax))))
  N["M", , ] <- matrix(census$male[index], nrow = length(years))
  N["F", , ] <- matrix(census$female[index], nrow = length(years))
  prob_repro <- rbind(M = spd[["male-prob-repro"]], F = spd[["fem-prob-repro"]])
  output <- rbind(M = spd[["male-asrp"]], F = spd[["fem-asrf"]])
  weights <- prob_repro * output
  E <- N
  reproductive_ages <- data.frame(sex = c("M", "F"), min_age = NA_integer_, max_age = NA_integer_)
  for (sex in c("M", "F")) {
    ages <- which(weights[sex, ] > 0)
    if (!length(ages)) stop("No reproductive ages for sex ", sex, ".", call. = FALSE)
    reproductive_ages[reproductive_ages$sex == sex, c("min_age", "max_age")] <- list(min(ages), max(ages))
    for (year in as.character(years)) {
      contribution <- N[sex, year, ] * weights[sex, ]
      total <- sum(contribution)
      if (!is.finite(total)) stop("Nonfinite reproductive contribution for sex ", sex, " in year ", year, ".", call. = FALSE)
      E[sex, year, ] <- if (total > 0) contribution / total else NA_real_
    }
  }
  list(N = N, E = E, years = years, pop = pop, reproductive_ages = reproductive_ages,
    prob_repro = prob_repro, relative_output = output,
    survival_into_age = rbind(M = spd[["male-surv-probs"]], F = spd[["fem-surv-probs"]]))
}
