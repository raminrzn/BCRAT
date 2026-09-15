test_that("model_run returns risk for the sample cohort", {
  out <- model_run(get_sample_input())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 5L)
  expect_true(all(c("risk", "risk_percent") %in% names(out)))
  expect_true(all(out$risk > 0 & out$risk < 1))
  expect_equal(out$risk_percent, round(100 * out$risk, 2))
})

test_that("model_run accepts a single-person list and the default", {
  out <- model_run(get_default_input())
  expect_equal(nrow(out), 1L)
  expect_true(out$risk > 0 && out$risk < 1)
  expect_equal(model_run(NULL)$risk, out$risk)  # NULL -> default
})

test_that("model_run ignores unknown fields but requires the core predictors", {
  ok <- get_default_input()
  ok$some_extra <- 1                     # extra field -> ignored, not fatal
  expect_no_error(model_run(ok))

  short <- get_default_input()
  short$age <- NULL
  expect_error(model_run(short), "Missing required variable")
})

test_that("model_run accepts unwrapped named args (do.call style)", {
  d <- get_default_input()
  expect_equal(model_run(d)$risk, do.call(model_run, d)$risk)
})

test_that("model_run accepts NCI-style raw field-name aliases", {
  viaAlias <- model_run(list(T1 = 50, race = "white", AgeMen = 13, Age1st = 24,
                             N_Biop = 0, N_Rels = 0, HypPlas = "unknown"))
  viaCanon <- model_run(list(age = 50, race = "white", age_menarche = 13,
                             age_first_birth = 24, n_biopsies = 0, n_relatives = 0,
                             atypical_hyperplasia = "unknown"))
  expect_equal(viaAlias$risk, viaCanon$risk)

  viaAlias2 <- model_run(list(age = 45, race = "hispanic", age_menarche = 12,
                              age_first_birth = 26, n_biopsies = 0, n_relatives = 1,
                              nativity = "foreign_born"))
  viaCanon2 <- model_run(list(age = 45, race = "hispanic", age_menarche = 12,
                              age_first_birth = 26, n_biopsies = 0, n_relatives = 1,
                              hispanic_nativity = "foreign_born"))
  expect_equal(viaAlias2$risk, viaCanon2$risk)
})

test_that("get_sample_input(n) limits rows and validates n", {
  expect_equal(nrow(get_sample_input(2)), 2L)
  expect_error(get_sample_input(0), "positive integer")
})

test_that("a JSON null is treated as an absent field, not a fatal input", {
  # hispanic_nativity/asian_subgroup are NA for anyone who isn't Hispanic or
  # Asian, so they serialise as null and come back as NULL on the next call.
  # Before this was handled, as.data.frame() died with "differing number of
  # rows" and OpenCPU surfaced it as an opaque HTTP 400.
  payload <- list(age = 50, race = "white", age_menarche = 13,
                  age_first_birth = 24, n_biopsies = 0, n_relatives = 0,
                  atypical_hyperplasia = "unknown",
                  hispanic_nativity = NULL, asian_subgroup = NULL)
  expect_no_error(out <- model_run(payload))
  expect_equal(out$risk, bcrat(50, "white", 13, 24, 0, 0), tolerance = 1e-12)
})

test_that("the model can re-consume its own serialised output", {
  # The platform builds its example payload from a previous response, so
  # gateway output must survive a round trip back through model_run().
  js <- gateway(func = "model_run", model_input = get_sample_input())
  round_tripped <- jsonlite::fromJSON(js, simplifyVector = FALSE)
  for (row in round_tripped) {
    row$risk <- NULL
    row$risk_percent <- NULL
    expect_no_error(model_run(row))
  }
})
