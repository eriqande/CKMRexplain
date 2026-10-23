library(CKMRexplain)

ages <- 1:3
spd <- list("max-age" = 3, "male-prob-repro" = c(1, 0, 1),
  "fem-prob-repro" = c(0, 1, 1), "male-asrp" = ages^2,
  "fem-asrf" = ages^2, "male-surv-probs" = rep(0.85, 3),
  "fem-surv-probs" = rep(0.85, 3))
census <- expand.grid(year = 0:20, age = ages, pop = 0)
census$male <- 10 + census$year + census$age
census$female <- 20 + 2 * census$year + census$age
inputs <- cohort_demography_from_spd(spd, census)
stopifnot(identical(dim(inputs$N), c(2L, 21L, 3L)), inputs$N["F", "10", "2"] == 42)
out <- pedigree_cohort_probabilities(inputs, 15, 17, 3)
probs <- out$cohort_probabilities
mass <- aggregate(probability ~ individual, probs, sum)
stopifnot(nrow(mass) == 30, all(abs(mass$probability - 1) < 1e-12))
father <- subset(probs, individual == "A-2")
stopifnot(father$probability[father$cohort == 13] == 0,
  isTRUE(all.equal(father$probability[father$cohort == 12], 28 * 9 / (26 + 28 * 9))))

# Enumerate all parental ages along each specified path independently.
for (index in 2:15) {
  path <- character()
  i <- index
  while (i > 1) {
    path <- c(if (i %% 2 == 0) "M" else "F", path)
    i <- floor(i / 2)
  }
  branches <- data.frame(cohort = 15, probability = 1)
  for (sex in path) {
    branches <- do.call(rbind, lapply(seq_len(nrow(branches)), function(j) {
      year <- branches$cohort[j]
      data.frame(cohort = year - ages,
        probability = branches$probability[j] * inputs$E[sex, as.character(year), ])
    }))
  }
  brute <- aggregate(probability ~ cohort, branches, sum)
  actual <- probs[probs$individual == paste0("A-", index), ]
  expected <- brute$probability[match(actual$cohort, brute$cohort)]
  stopifnot(max(abs(actual$probability - expected)) < 1e-12)
}
scaled <- census
scaled$male <- scaled$male / (scaled$year + 1)
scaled$female <- scaled$female / (scaled$year + 2)
stopifnot(isTRUE(all.equal(inputs$E, cohort_demography_from_spd(spd, scaled)$E)))
survival_changed <- spd
survival_changed$`male-surv-probs` <- rep(1, 3)
stopifnot(identical(inputs$E, cohort_demography_from_spd(survival_changed, census)$E))

error_message <- function(expr) tryCatch({ force(expr); stop("Expected an error") }, error = conditionMessage)
stopifnot(grepl("minimum cohort is 9", error_message(pedigree_cohort_probabilities(inputs, 8, 17, 3))),
  grepl("after the latest census", error_message(pedigree_cohort_probabilities(inputs, 15, 17, 7))),
  grepl("supported sample-cohort range", error_message(pedigree_cohort_probabilities(inputs, 21, 17, 3))),
  grepl("Incomplete census", error_message(cohort_demography_from_spd(spd, census[-1, ]))),
  grepl("Duplicate census", error_message(cohort_demography_from_spd(spd, rbind(census, census[1, ])))))
zero <- census
zero$male[zero$year == 15] <- 0
stopifnot(grepl("Undefined reproductive distribution", error_message(
  pedigree_cohort_probabilities(cohort_demography_from_spd(spd, zero), 15, 17, 3))))
stopifnot(nrow(pedigree_cohort_probabilities(inputs, 0, 20, 0)$cohort_probabilities) == 2)
