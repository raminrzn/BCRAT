# Reference values were NOT hand-derived: they were produced by independently
# executing NCI's own reference algorithm -- the CRAN "BCRA" package (v2.1.2,
# GPL >= 2), which implements the identical published BCRAT algorithm --
# either from its own bundled cross-validation fixtures (via the independent
# epiverse/bcra-js port's R-generated
# test/fixtures/r-reference/race-specific.json, itself produced by running
# the CRAN package) or by running the CRAN package's source directly for the
# additional edge cases below. All rows use a fixed 5-year window (T2 = T1+5).

# Three raw input profiles shared across all eleven NCI race/ethnicity codes.
# Ages 38/40/42 (all within this package's supported 35-85 range); p1/p3 were
# re-derived at these ages from the reference source directly (the original
# fixture file used ages 25/30/40, which exercise the underlying algorithm's
# wider 20-90 mathematical domain but fall outside the 35-85 range this
# package enforces, matching NCI's own deployed consumer tool).
profiles <- list(
  p1 = list(age = 38, n_biopsies = 0, atypical_hyperplasia = "unknown",
            age_menarche = 14, age_first_birth = 22, n_relatives = 0),
  p2 = list(age = 40, n_biopsies = 1, atypical_hyperplasia = "no",
            age_menarche = 11, age_first_birth = "nulliparous", n_relatives = 1),
  p3 = list(age = 42, n_biopsies = 0, atypical_hyperplasia = "unknown",
            age_menarche = 13, age_first_birth = "unknown", n_relatives = 0)
)

groups <- data.frame(
  race     = c("white", "black", "hispanic", "american_indian_or_alaska_native", "hispanic",
               "asian", "asian", "asian", "asian", "asian", "asian"),
  nativity = c(NA, NA, "us_born", NA, "foreign_born", NA, NA, NA, NA, NA, NA),
  subgroup = c(NA, NA, NA, NA, NA, "chinese", "japanese", "filipino", "hawaiian",
               "other_pacific_islander", "other_asian"),
  p1 = c(0.3668333224, 0.3568809234, 0.4592488062, 0.3668333224, 0.2050263136,
         0.2329454751, 0.2840942254, 0.2803548156, 0.4176570317, 0.2162906768, 0.2052195495),
  p2 = c(1.8983595395, 1.0088314365, 0.8746778404, 1.8983595395, 2.5618113231,
         1.5498436810, 1.9500749503, 1.9115693839, 2.8010693030, 1.2788445091, 1.3434769118),
  p3 = c(0.4755871335, 0.6897035991, 0.5560259079, 0.4755871335, 0.2550309452,
         0.2902748281, 0.3661401045, 0.3579710683, 0.5133473567, 0.1933365983, 0.2283614340),
  stringsAsFactors = FALSE
)

test_that("bcrat reproduces the NCI reference algorithm across all race/ethnicity groups", {
  for (g in seq_len(nrow(groups))) {
    for (p in c("p1", "p2", "p3")) {
      prof <- profiles[[p]]
      got <- bcrat(
        age = prof$age, race = groups$race[g],
        age_menarche = prof$age_menarche, age_first_birth = prof$age_first_birth,
        n_biopsies = prof$n_biopsies, n_relatives = prof$n_relatives,
        atypical_hyperplasia = prof$atypical_hyperplasia,
        hispanic_nativity = groups$nativity[g], asian_subgroup = groups$subgroup[g]
      )
      label <- paste(groups$race[g], groups$nativity[g], groups$subgroup[g], p)
      expect_equal(100 * got, groups[[p]][g], tolerance = 1e-6, label = label)
    }
  }
})

test_that("American Indian/Alaska Native matches White exactly (NCI substitutes the White model)", {
  args <- list(age_menarche = 13, age_first_birth = 24, n_biopsies = 1,
               n_relatives = 0, atypical_hyperplasia = "no")
  white <- do.call(bcrat, c(list(age = 48, race = "white"), args))
  aian  <- do.call(bcrat, c(list(age = 48, race = "american_indian_or_alaska_native"), args))
  expect_equal(aian, white)
})

test_that("boundary cases: straddling age 50, >=50-only, and the valid age extremes", {
  # White, straddling age 50 (ages 48-52 inclusive contribute to the 5-yr sum)
  expect_equal(100 * bcrat(48, "white", 13, 24, 1, 0, "no"), 1.0745946093, tolerance = 1e-6)
  # White, entirely >=50, with atypical hyperplasia
  expect_equal(100 * bcrat(50, "white", 11, 22, 2, 1, "yes"), 5.8684594548, tolerance = 1e-6)
  # Black, entirely >=50
  expect_equal(100 * bcrat(55, "black", 10, "unknown", 2, 2, "yes"), 7.2132985085, tolerance = 1e-6)
  # Foreign-born Hispanic, straddling age 50
  expect_equal(100 * bcrat(47, "hispanic", 12, 26, 0, 1, "unknown", hispanic_nativity = "foreign_born"),
               1.4007685831, tolerance = 1e-6)
  # Chinese, entirely >=50
  expect_equal(100 * bcrat(55, "asian", 14, 25, 1, 0, "unknown", asian_subgroup = "chinese"),
               1.4787809518, tolerance = 1e-6)
  # American Indian/Alaska Native, straddling age 50 -- must equal the White figure above
  expect_equal(100 * bcrat(48, "american_indian_or_alaska_native", 13, 24, 1, 0, "no"),
               1.0745946093, tolerance = 1e-6)
  # White at the minimum valid age (35)
  expect_equal(100 * bcrat(35, "white", 13, 25, 0, 0, "unknown"), 0.3242306364, tolerance = 1e-6)
  # White at the maximum valid age (85)
  expect_equal(100 * bcrat(85, "white", 13, 25, 1, 0, "no"), 1.6857511273, tolerance = 1e-6)
})

test_that("Hispanic-specific pooling (2+ biopsies -> 1) vs. AIAN/White (no pooling)", {
  # Same raw inputs (4 biopsies, nulliparous/25yo first birth), age 35-40.
  expect_equal(100 * bcrat(35, "hispanic", 11, 25, 4, 0, hispanic_nativity = "us_born"),
               0.3162827, tolerance = 1e-6)
  expect_equal(100 * bcrat(35, "american_indian_or_alaska_native", 11, "nulliparous", 4, 0),
               1.0229802, tolerance = 1e-6)
})

test_that("age_first_birth is inert for Black (CARE model excludes it)", {
  base <- bcrat(45, "black", 13, 20, 1, 1, "no")
  alt  <- bcrat(45, "black", 13, 35, 1, 1, "no")
  expect_equal(base, alt)
})

test_that("age_menarche is inert for US-born Hispanic (coefficient is 0)", {
  base <- bcrat(45, "hispanic", 11, 24, 1, 1, "no", hispanic_nativity = "us_born")
  alt  <- bcrat(45, "hispanic", 14, 24, 1, 1, "no", hispanic_nativity = "us_born")
  expect_equal(base, alt)
})

test_that("nulliparous and unknown sentinels are accepted for age_first_birth/age_menarche/etc.", {
  expect_no_error(bcrat(50, "white", "unknown", "nulliparous", "unknown", "unknown"))
  expect_no_error(bcrat(50, "white", 13, 24, 0, 0, atypical_hyperplasia = "unknown"))
})

test_that("numeric race codes 0-5 match their character equivalents", {
  args <- list(age_menarche = 13, age_first_birth = 24, n_biopsies = 0, n_relatives = 0)
  expect_equal(do.call(bcrat, c(list(age = 50, race = 0), args)),
               do.call(bcrat, c(list(age = 50, race = "white"), args)))
  expect_equal(do.call(bcrat, c(list(age = 50, race = 1), args)),
               do.call(bcrat, c(list(age = 50, race = "white"), args)))
  expect_equal(do.call(bcrat, c(list(age = 50, race = 2), args)),
               do.call(bcrat, c(list(age = 50, race = "black"), args)))
  expect_equal(do.call(bcrat, c(list(age = 50, race = 5), args)),
               do.call(bcrat, c(list(age = 50, race = "american_indian_or_alaska_native"), args)))
})

test_that("a mixed cohort (not all Hispanic/Asian) vectorises without spurious errors", {
  # Regression test: hispanic_nativity/asian_subgroup carry NA on rows where
  # they don't apply; resolving them must not validate those irrelevant NAs.
  out <- bcrat(
    age = c(50, 45, 62, 58, 70),
    race = c("white", "black", "hispanic", "asian", "american_indian_or_alaska_native"),
    age_menarche = c(13, 12, 14, 11, 13),
    age_first_birth = c(24, "nulliparous", 26, 22, 20),
    n_biopsies = c(0, 1, 0, 2, 1),
    n_relatives = c(0, 1, 1, 0, 2),
    atypical_hyperplasia = c("unknown", "no", "unknown", "yes", "no"),
    hispanic_nativity = c(NA, NA, "foreign_born", NA, NA),
    asian_subgroup = c(NA, NA, NA, "chinese", NA)
  )
  expect_length(out, 5L)
  expect_true(all(out > 0 & out < 1))
  # Row 3 (Hispanic, foreign-born) must match the equivalent single-row call.
  expect_equal(out[3], bcrat(62, "hispanic", 14, 26, 0, 1, "unknown", hispanic_nativity = "foreign_born"))
  # Row 4 (Asian, Chinese) must match the equivalent single-row call.
  expect_equal(out[4], bcrat(58, "asian", 11, 22, 2, 0, "yes", asian_subgroup = "chinese"))
})

test_that("invalid inputs are rejected with informative errors", {
  base <- list(age_menarche = 13, age_first_birth = 24, n_biopsies = 0, n_relatives = 0)
  expect_error(do.call(bcrat, c(list(age = 34, race = "white"), base)), "35 to 85")
  expect_error(do.call(bcrat, c(list(age = 86, race = "white"), base)), "35 to 85")
  expect_error(do.call(bcrat, c(list(age = 50.5, race = "white"), base)), "35 to 85")
  expect_error(do.call(bcrat, c(list(age = 50, race = "martian"), base)), "Unrecognised .race.")
  expect_error(do.call(bcrat, c(list(age = 50, race = 9), base)), "0-5")
  expect_error(
    bcrat(50, "hispanic", 13, 24, 0, 0, hispanic_nativity = "unsure"),
    "Unrecognised .hispanic_nativity."
  )
  expect_error(
    bcrat(50, "asian", 13, 24, 0, 0, asian_subgroup = "vietnamese"),
    "Unrecognised .asian_subgroup."
  )
  expect_error(do.call(bcrat, c(list(age = 50, race = "white"), list(
    age_menarche = 13, age_first_birth = 24, n_biopsies = -1, n_relatives = 0))),
    "n_biopsies")
})

test_that("bcrat is vectorised and recycles length-1 arguments", {
  out <- bcrat(age = c(50, 60), race = "white", age_menarche = 13,
               age_first_birth = 24, n_biopsies = 0, n_relatives = 0)
  expect_length(out, 2L)
  expect_true(out[2] > out[1])  # risk rises with age, all else equal
})
