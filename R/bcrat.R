# ---------------------------------------------------------------------------
# BCRAT — NCI Breast Cancer Risk Assessment Tool ("Gail model") — core
#
# Gail MH, Brinton LA, Byar DP, et al. Projecting individualized probabilities
#   of developing breast cancer for white females who are being examined
#   annually. J Natl Cancer Inst. 1989;81(24):1879-1886. doi:10.1093/jnci/81.24.1879
# Costantino JP, Gail MH, Pee D, et al. Validation studies for models
#   projecting the risk of invasive and total breast cancer incidence.
#   J Natl Cancer Inst. 1999;91(18):1541-1548. doi:10.1093/jnci/91.18.1541
# Gail MH, Costantino JP, Pee D, et al. Projecting individualized absolute
#   invasive breast cancer risk in African American women. J Natl Cancer Inst.
#   2007;99(23):1782-1792. doi:10.1093/jnci/djm223
# Matsuno RK, Costantino JP, Ziegler RG, et al. Projecting individualized
#   absolute invasive breast cancer risk in Asian and Pacific Islander American
#   women. J Natl Cancer Inst. 2011;103(12):951-961. doi:10.1093/jnci/djr154
# Banegas MP, John EM, Slattery ML, et al. Projecting individualized absolute
#   invasive breast cancer risk in US Hispanic women. J Natl Cancer Inst.
#   2017;109(2):djw215. doi:10.1093/jnci/djw215
#
# The model combines a race-specific relative-risk (logistic) sub-model with
# race-specific baseline hazards (SEER invasive-breast-cancer incidence and
# NCHS competing, non-breast-cancer mortality) in a discrete, single-year
# competing-risk life table (Gail et al. 1989 eq. 1; Benichou & Gail,
# Biometrics 1990;46:813-826). Coefficients and rate tables below were
# cross-verified by independently executing NCI's own reference algorithm,
# distributed under GPL (>= 2) as the CRAN package "BCRA"
# (https://cran.r-project.org/package=BCRA, maintained by NCI's own BCRAT
# technical-support contact), against published worked examples in the papers
# above and against WhistlerBrown/epiverse's independent "bcra-js" port
# (GPL-3.0-or-later). This implementation is written from scratch against the
# published algorithm and coefficients; it does not reuse BCRA's source code.
# ---------------------------------------------------------------------------

# Relative-risk coefficients (log-RR scale): n_biopsies, age_menarche,
# age_first_birth, n_relatives, age>=50 x n_biopsies, age_first_birth x
# n_relatives. All six Asian/Pacific Islander subgroups share one coefficient
# vector (relative risks were tested for homogeneity across subgroups in
# Matsuno et al. 2011); American Indian/Alaska Native has no independently
# validated model and uses the White coefficients (see .bcrat_model_group()).
.bcrat_beta <- rbind(
  white                 = c(0.5292641686,     0.0940103059,     0.2186262218,     0.9583027845,     -0.2880424830, -0.1908113865),
  black                 = c(0.1822121131,     0.2672530336,     0.0,              0.4757242578,     -0.1119411682,  0.0),
  hispanic_us_born      = c(0.0970783641,     0.0,              0.2318368334,     0.1666854410,      0.0,           0.0),
  hispanic_foreign_born = c(0.4798624017,     0.2593922322,     0.4669246218,     0.9076679727,      0.0,           0.0),
  asian                 = c(0.55263612260619, 0.07499257592975, 0.27638268294593, 0.79185633720481,  0.0,           0.0)
)
colnames(.bcrat_beta) <- c("n_biopsies", "age_menarche", "age_first_birth",
                           "n_relatives", "age50_x_biopsies", "first_birth_x_relatives")

# Baseline invasive breast cancer incidence (lambda1) and competing,
# non-breast-cancer mortality (lambda2): 14 five-year age bands, [20,25) ...
# [85,90). White = SEER 1983-87 / NCHS 1985-87; Black = SEER 1994-98 / NCHS
# 1996-2000; Hispanic = California Cancer Registry/SEER 1995-2004 (US-born and
# foreign-born are separately tabulated); Asian subgroups = SEER-18 / NCHS
# 1998-2002. American Indian/Alaska Native uses the White tables (see below).
.bcrat_lambda1 <- rbind(
  white                  = c(0.00001000, 0.00007600, 0.00026600, 0.00066100, 0.00126500, 0.00186600, 0.00221100, 0.00272100, 0.00334800, 0.00392300, 0.00417800, 0.00443900, 0.00442100, 0.00410900),
  black                  = c(0.00002696, 0.00011295, 0.00031094, 0.00067639, 0.00119444, 0.00187394, 0.00241504, 0.00291112, 0.00310127, 0.00366560, 0.00393132, 0.00408951, 0.00396793, 0.00363712),
  hispanic_us_born       = c(0.0000166,  0.0000741,  0.0002740,  0.0006099,  0.0012225,  0.0019027,  0.0023142,  0.0028357,  0.0031144,  0.0030794,  0.0033344,  0.0035082,  0.0025308,  0.0020414),
  hispanic_foreign_born  = c(0.0000102,  0.0000531,  0.0001578,  0.0003602,  0.0007617,  0.0011599,  0.0014111,  0.0017245,  0.0020619,  0.0023603,  0.0025575,  0.0028227,  0.0028295,  0.0025868),
  chinese                = c(0.000004059636, 0.000045944465, 0.000188279352, 0.000492930493, 0.000913603501, 0.001471537353, 0.001421275482, 0.001970946494, 0.001674745804, 0.001821581075, 0.001834477198, 0.001919911972, 0.002233371071, 0.002247315779),
  japanese               = c(0.000000000001, 0.000099483924, 0.000287041681, 0.000545285759, 0.001152211095, 0.001859245108, 0.002606291272, 0.003221751682, 0.004006961859, 0.003521715275, 0.003593038294, 0.003589303081, 0.003538507159, 0.002051572909),
  filipino               = c(0.000007500161, 0.000081073945, 0.000227492565, 0.000549786433, 0.001129400541, 0.001813873795, 0.002223665639, 0.002680309266, 0.002891219230, 0.002534421279, 0.002457159409, 0.002286616920, 0.001814802825, 0.001750879130),
  hawaiian               = c(0.000045080582, 0.000098570724, 0.000339970860, 0.000852591429, 0.001668562761, 0.002552703284, 0.003321774046, 0.005373001776, 0.005237808549, 0.005581732512, 0.005677419355, 0.006513409962, 0.003889457523, 0.002949061662),
  other_pacific_islander = c(0.000000000001, 0.000071525212, 0.000288799028, 0.000602250698, 0.000755579402, 0.000766406354, 0.001893124938, 0.002365580107, 0.002843933070, 0.002920921732, 0.002330395655, 0.002036291235, 0.001482683983, 0.001012248203),
  other_asian            = c(0.000012355409, 0.000059526456, 0.000184320831, 0.000454677273, 0.000791265338, 0.001048462801, 0.001372467817, 0.001495473711, 0.001646746198, 0.001478363563, 0.001216010125, 0.001067663700, 0.001376104012, 0.000661576644)
)

.bcrat_lambda2 <- rbind(
  white                  = c(0.00049300, 0.00053100, 0.00062500, 0.00082500, 0.00130700, 0.00218100, 0.00365500, 0.00585200, 0.00943900, 0.01502800, 0.02383900, 0.03883200, 0.06682800, 0.14490800),
  black                  = c(0.00074354, 0.00101698, 0.00145937, 0.00215933, 0.00315077, 0.00448779, 0.00632281, 0.00963037, 0.01471818, 0.02116304, 0.03266035, 0.04564087, 0.06835185, 0.13271262),
  hispanic_us_born       = c(0.0003561,  0.0004038,  0.0005281,  0.0008875,  0.0013987,  0.0020769,  0.0030912,  0.0046960,  0.0076050,  0.0120555,  0.0193805,  0.0288386,  0.0429634,  0.0740349),
  hispanic_foreign_born  = c(0.0003129,  0.0002908,  0.0003515,  0.0004943,  0.0007807,  0.0012840,  0.0020325,  0.0034533,  0.0058674,  0.0096888,  0.0154429,  0.0254675,  0.0448037,  0.1125678),
  chinese                = c(0.000210649076, 0.000192644865, 0.000244435215, 0.000317895949, 0.000473261994, 0.000800271380, 0.001217480226, 0.002099836508, 0.003436889186, 0.006097405623, 0.010664526765, 0.020148678452, 0.037990796590, 0.098333900733),
  japanese               = c(0.000173593803, 0.000295805882, 0.000228322534, 0.000363242389, 0.000590633044, 0.001086079485, 0.001859999966, 0.003216600974, 0.004719402141, 0.008535331402, 0.012433511681, 0.020230197885, 0.037725498348, 0.106149118663),
  filipino               = c(0.000229120979, 0.000262988494, 0.000314844090, 0.000394471908, 0.000647622610, 0.001170202327, 0.001809380379, 0.002614170568, 0.004483330681, 0.007393665092, 0.012233059675, 0.021127058106, 0.037936954809, 0.085138518334),
  hawaiian               = c(0.000563507269, 0.000369640217, 0.001019912579, 0.001234013911, 0.002098344078, 0.002982934175, 0.005402445702, 0.009591474245, 0.016315472607, 0.020152229069, 0.027354838710, 0.050446998723, 0.072262026612, 0.145844504021),
  other_pacific_islander = c(0.000465500812, 0.000600466920, 0.000851057138, 0.001478265376, 0.001931486788, 0.003866623959, 0.004924932309, 0.008177071806, 0.008638202890, 0.018974658371, 0.029257567105, 0.038408980974, 0.052869579345, 0.074745721133),
  other_asian            = c(0.000212632332, 0.000242170741, 0.000301552711, 0.000369053354, 0.000543002943, 0.000893862331, 0.001515172239, 0.002574669551, 0.004324370426, 0.007419621918, 0.013251765130, 0.022291427490, 0.041746550635, 0.087485802065)
)

# 1 minus the attributable-risk fraction, [age<50, age>=50]. Multiplies the
# incidence hazard alongside the relative risk (Gail et al. 1989 eq. 1).
.bcrat_far <- rbind(
  white                 = c(0.5788413,      0.5788413),
  black                 = c(0.72949880,     0.74397137),
  hispanic_us_born      = c(0.749294788397, 0.778215491668),
  hispanic_foreign_born = c(0.428864989813, 0.450352338746),
  asian                 = c(0.47519806426735, 0.50316401683903)
)

# American Indian/Alaska Native (and, in NCI's own tool, "unknown" race) has no
# independently validated relative-risk or baseline-rate model; the deployed
# tool substitutes the White model in full. Asian/Pacific Islander subgroups
# share one relative-risk/attributable-risk model but have distinct baseline
# rate tables (.bcrat_lambda1/2 above already carry six separate rows).
.bcrat_model_group <- function(group) {
  ifelse(group == "american_indian_or_alaska_native", "white",
  ifelse(group %in% c("chinese", "japanese", "filipino", "hawaiian",
                       "other_pacific_islander", "other_asian"),
         "asian", group))
}
.bcrat_lambda_group <- function(group) {
  ifelse(group == "american_indian_or_alaska_native", "white", group)
}

# Five-year age band index (1 = [20,25), ..., 14 = [85,90)), clamped.
.bcrat_band <- function(age) pmin(pmax(floor((age - 20) / 5) + 1L, 1L), 14L)

# Race/ethnicity aliases -> canonical category. "Unknown" is mapped the same
# way NCI's own tool maps it: to the American Indian/Alaska Native ("Other")
# bucket, i.e. the White model.
.bcrat_race_aliases <- c(
  "white" = "white", "caucasian" = "white", "non-hispanic white" = "white",
  "black" = "black", "african american" = "black", "african-american" = "black",
  "hispanic" = "hispanic", "latino" = "hispanic", "latina" = "hispanic",
  "asian" = "asian", "asian american" = "asian", "asian-american" = "asian",
  "asian/pacific islander" = "asian", "asian or pacific islander" = "asian",
  "american indian or alaska native" = "american_indian_or_alaska_native",
  "american indian or alaskan native" = "american_indian_or_alaska_native",
  "american_indian_or_alaska_native" = "american_indian_or_alaska_native",
  "american indian" = "american_indian_or_alaska_native",
  "alaska native" = "american_indian_or_alaska_native",
  "alaskan native" = "american_indian_or_alaska_native",
  "indigenous" = "american_indian_or_alaska_native",
  "unknown" = "american_indian_or_alaska_native"
)
.bcrat_race_codes <- c("white", "black", "hispanic", "asian", "american_indian_or_alaska_native")

.bcrat_nativity_aliases <- c(
  "us_born" = "us_born", "us-born" = "us_born", "us born" = "us_born",
  "us" = "us_born", "united states" = "us_born", "in the us" = "us_born",
  "foreign_born" = "foreign_born", "foreign-born" = "foreign_born",
  "foreign born" = "foreign_born", "foreign" = "foreign_born",
  "outside us" = "foreign_born", "outside the us" = "foreign_born"
)

.bcrat_subgroup_aliases <- c(
  "chinese" = "chinese",
  "japanese" = "japanese",
  "filipino" = "filipino",
  "hawaiian" = "hawaiian", "native hawaiian" = "hawaiian",
  "other_pacific_islander" = "other_pacific_islander",
  "other pacific islander" = "other_pacific_islander",
  "pacific islander" = "other_pacific_islander",
  "other_asian" = "other_asian", "other asian" = "other_asian"
)

.bcrat_lookup_alias <- function(x, alias_map, arg_name, choices_msg) {
  key <- tolower(trimws(as.character(x)))
  canon <- unname(alias_map[key])
  if (any(is.na(canon))) {
    bad <- unique(x[is.na(canon)])
    stop("Unrecognised `", arg_name, "` value(s): ", paste(bad, collapse = ", "),
         ". ", choices_msg, call. = FALSE)
  }
  canon
}

.bcrat_subgroup_choices_msg <- paste(
  'Use one of "chinese", "japanese", "filipino", "hawaiian",',
  '"other_pacific_islander", "other_asian".'
)

# Resolve `race` (character alias or numeric 0-5) to a canonical top-level
# category. 0 and 1 both mean White (0 accepted as a synonym, matching the
# other ModelsCloud packages in this family).
.bcrat_resolve_race <- function(race) {
  if (is.numeric(race)) {
    race[race == 0] <- 1L
    if (any(!race %in% seq_along(.bcrat_race_codes), na.rm = TRUE)) {
      stop("Numeric `race` codes must be 0-5: 0/1 White, 2 Black, 3 Hispanic, ",
           "4 Asian, 5 American Indian or Alaska Native.", call. = FALSE)
    }
    return(.bcrat_race_codes[race])
  }
  .bcrat_lookup_alias(race, .bcrat_race_aliases, "race", "See ?bcrat for accepted categories.")
}

# Resolve `race` (+ `hispanic_nativity` / `asian_subgroup` where relevant) to
# one of the eleven underlying model groups. Nativity/subgroup are validated
# only on the rows where they actually apply, so an NA/absent placeholder on
# an irrelevant row (e.g. a mixed cohort where most women aren't Hispanic)
# never triggers a spurious "unrecognised value" error.
.bcrat_resolve_group <- function(race, hispanic_nativity, asian_subgroup) {
  canon_race <- .bcrat_resolve_race(race)
  group <- canon_race

  is_hisp  <- canon_race == "hispanic"
  is_asian <- canon_race == "asian"

  if (any(is_hisp)) {
    nat <- if (is.null(hispanic_nativity)) rep("us_born", sum(is_hisp)) else hispanic_nativity[is_hisp]
    nat[is.na(nat)] <- "us_born"
    nat <- .bcrat_lookup_alias(nat, .bcrat_nativity_aliases, "hispanic_nativity",
                               'Use "us_born" or "foreign_born".')
    group[is_hisp] <- ifelse(nat == "us_born", "hispanic_us_born", "hispanic_foreign_born")
  }
  if (any(is_asian)) {
    sub <- if (is.null(asian_subgroup)) rep("other_asian", sum(is_asian)) else asian_subgroup[is_asian]
    sub[is.na(sub)] <- "other_asian"
    sub <- .bcrat_lookup_alias(sub, .bcrat_subgroup_aliases, "asian_subgroup",
                               .bcrat_subgroup_choices_msg)
    group[is_asian] <- sub
  }
  group
}

# Number of prior breast biopsies -> category 0/1/2+ (99 or "unknown" -> 0,
# treated as no biopsies). Hispanic (either nativity) pools 2+ with 1 (San
# Francisco Bay Area Breast Cancer Study recoding, Banegas et al. 2017).
.bcrat_recode_nb <- function(n_biopsies, group) {
  is_unk <- is.character(n_biopsies) & tolower(trimws(n_biopsies)) %in% c("unknown", "unk", "99")
  x <- suppressWarnings(as.numeric(ifelse(is_unk, NA, n_biopsies)))
  if (any(is.na(x) & !is_unk)) {
    stop('`n_biopsies` must be a non-negative number or "unknown".', call. = FALSE)
  }
  if (any(x < 0, na.rm = TRUE)) {
    stop('`n_biopsies` must be a non-negative number or "unknown".', call. = FALSE)
  }
  cat <- ifelse(is_unk | x == 0, 0L, ifelse(x == 1, 1L, 2L))
  is_hisp <- group %in% c("hispanic_us_born", "hispanic_foreign_born")
  cat[is_hisp & cat == 2L] <- 1L
  cat
}

# Age at menarche -> category 0 (>=14 or unknown), 1 (12-13), 2 (<12). Black
# collapses category 2 into 1 (Gail et al. 2007). (US-born Hispanic's
# coefficient for this term is 0, so its category value never affects risk.)
.bcrat_recode_am <- function(age_menarche, group) {
  is_unk <- is.character(age_menarche) & tolower(trimws(age_menarche)) %in% c("unknown", "unk", "99")
  x <- suppressWarnings(as.numeric(ifelse(is_unk, NA, age_menarche)))
  if (any(is.na(x) & !is_unk)) {
    stop('`age_menarche` must be an age in years or "unknown".', call. = FALSE)
  }
  cat <- ifelse(is_unk | x >= 14, 0L, ifelse(x >= 12, 1L, 2L))
  cat[group == "black" & cat == 2L] <- 1L
  cat
}

# Age at first live birth -> category 0 (<20 or unknown), 1 (20-24), 2 (25-29
# or nulliparous), 3 (>=30). Hispanic (either nativity) demotes true 25-29
# into category 1 (merged with 20-24) and >=30 into category 2 (merged with
# nulliparous) (Banegas et al. 2017). (Black's coefficient and its interaction
# with n_relatives are both 0, so age_first_birth never affects Black risk —
# consistent with the CARE model excluding this covariate entirely.)
.bcrat_recode_af <- function(age_first_birth, group) {
  is_null <- is.character(age_first_birth) &
    tolower(trimws(age_first_birth)) %in% c("nulliparous", "no births", "none", "nullip", "98")
  is_unk <- is.character(age_first_birth) &
    tolower(trimws(age_first_birth)) %in% c("unknown", "unk", "99")
  x <- suppressWarnings(as.numeric(ifelse(is_null | is_unk, NA, age_first_birth)))
  if (any(is.na(x) & !is_null & !is_unk)) {
    stop('`age_first_birth` must be an age in years, "nulliparous", or "unknown".', call. = FALSE)
  }
  base <- ifelse(is_unk, 0L,
          ifelse(is_null, 2L,
          ifelse(x < 20, 0L,
          ifelse(x < 25, 1L,
          ifelse(x < 30, 2L, 3L)))))
  is_hisp <- group %in% c("hispanic_us_born", "hispanic_foreign_born")
  demote_25_29 <- is_hisp & base == 2L & !is_null
  demote_30p   <- is_hisp & base == 3L
  base[demote_25_29] <- 1L
  base[demote_30p]   <- 2L
  base
}

# Number of first-degree relatives with breast cancer -> category 0/1/2+ (99
# or "unknown" -> 0). Hispanic and all Asian/Pacific Islander subgroups pool
# 2+ with 1.
.bcrat_recode_nr <- function(n_relatives, group) {
  is_unk <- is.character(n_relatives) & tolower(trimws(n_relatives)) %in% c("unknown", "unk", "99")
  x <- suppressWarnings(as.numeric(ifelse(is_unk, NA, n_relatives)))
  if (any(is.na(x) & !is_unk)) {
    stop('`n_relatives` must be a non-negative number or "unknown".', call. = FALSE)
  }
  if (any(x < 0, na.rm = TRUE)) {
    stop('`n_relatives` must be a non-negative number or "unknown".', call. = FALSE)
  }
  cat <- ifelse(is_unk | x == 0, 0L, ifelse(x == 1, 1L, 2L))
  pools <- group %in% c("hispanic_us_born", "hispanic_foreign_born", "chinese",
                        "japanese", "filipino", "hawaiian",
                        "other_pacific_islander", "other_asian")
  cat[pools & cat == 2L] <- 1L
  cat
}

# Atypical-hyperplasia multiplier on the biopsy relative risk. Not applicable
# (multiplier 1.00) whenever the (possibly pooled) biopsy category is 0,
# regardless of what `atypical_hyperplasia` was passed.
.bcrat_r_hyp <- function(nb_cat, atypical_hyperplasia) {
  key <- tolower(trimws(as.character(atypical_hyperplasia)))
  hyp <- ifelse(key %in% c("yes", "y", "1", "true"), 1L,
         ifelse(key %in% c("no", "n", "0", "false"), 0L,
         ifelse(key %in% c("unknown", "unk", "na", "99"), 99L, NA_integer_)))
  if (any(is.na(hyp))) {
    stop('`atypical_hyperplasia` must be "yes", "no", or "unknown".', call. = FALSE)
  }
  ifelse(nb_cat == 0L, 1.00, ifelse(hyp == 0L, 0.93, ifelse(hyp == 1L, 1.82, 1.00)))
}

#' BCRAT 5-year invasive breast cancer risk
#'
#' Computes the 5-year probability of invasive breast cancer using the NCI
#' Breast Cancer Risk Assessment Tool (BCRAT), commonly known as the "Gail
#' model" in its Costantino et al. (1999) form, with the race/ethnicity-
#' specific extensions for Black (Gail et al. 2007), Hispanic (Banegas et al.
#' 2017), and Asian/Pacific Islander (Matsuno et al. 2011) women. All
#' arguments are vectorised and recycled to a common length.
#'
#' @details
#' The model combines a race-specific relative-risk sub-model with race-
#' specific baseline hazards in a single-year discrete competing-risk life
#' table (Gail et al. 1989, equation 1). For each of the five single years of
#' age from `age` to `age + 4`:
#'
#' \deqn{h = F(t) \times RR \times \lambda_1, \quad
#'       \tau = h + \lambda_2, \quad
#'       \pi = \frac{h}{\tau}\,e^{-\Sigma\tau}\left(1 - e^{-\tau}\right)}
#'
#' where \eqn{\lambda_1} and \eqn{\lambda_2} are the race-specific baseline
#' invasive-breast-cancer incidence and competing (non-breast-cancer)
#' mortality rates for the five-year age band containing that year, \eqn{RR}
#' is `RR1` (relative risk for ages <50) or `RR2` (ages >=50) as appropriate,
#' \eqn{F(t)} is one minus the attributable-risk fraction, and
#' \eqn{\Sigma\tau} is the running total hazard from all earlier years. The
#' 5-year risk is the sum of \eqn{\pi} over the five years.
#'
#' The relative risk is `RR1 = exp(LP1)` (ages <50) and
#' `RR2 = exp(LP1 + biopsies x age50_x_biopsies_beta)` (ages >=50), where
#'
#' \deqn{LP1 = \beta_{biop}\,NB + \beta_{men}\,AM + \beta_{birth}\,AF +
#'             \beta_{rel}\,NR + \beta_{birth \times rel}\,AF{\cdot}NR +
#'             \ln(R_{hyp})}
#'
#' `NB`, `AM`, `AF`, `NR` are the categorised biopsy, menarche, first-birth,
#' and relative-count covariates (see the parameter descriptions below for
#' exact category boundaries and race-specific pooling), and `R_hyp` is a
#' fixed multiplier for atypical hyperplasia among biopsied women (`0.93` no,
#' `1.82` yes, `1.00` unknown or not applicable).
#'
#' **American Indian/Alaska Native** (and "unknown" race) has no
#' independently validated relative-risk or baseline-rate model; NCI's
#' deployed tool substitutes the White model in full, and this package does
#' the same for exact parity.
#'
#' **Hispanic** requires a nativity (`hispanic_nativity`): US-born and
#' foreign-born Hispanic women have two entirely distinct fitted models (not
#' one model with a nativity adjustment), with materially different
#' coefficients and baseline rates. When not supplied, `"us_born"` is used.
#'
#' **Asian/Pacific Islander** requires a subgroup (`asian_subgroup`) for its
#' baseline rates (the six subgroups share one relative-risk model but have
#' distinct SEER incidence/mortality tables, and NCI's tool has no pooled
#' "combined" option). When not supplied, `"other_asian"` — the SEER category
#' for Asian ethnicities not individually tabulated — is used as the closest
#' available default; specify a subgroup for a more precise estimate.
#'
#' @param age Age in years, as a whole number from 35 to 85 (the range NCI's
#'   own tool supports; the underlying model is undefined outside it).
#' @param race Race/ethnicity, as a character label or numeric code.
#'   Accepted labels (case-insensitive): `"White"`, `"Black"`, `"Hispanic"`,
#'   `"Asian"`, `"American Indian or Alaska Native"` (plus common aliases such
#'   as `"Caucasian"`, `"African American"`, `"Latino"`, `"Unknown"`).
#'   Numeric codes: `1` (or `0`) White, `2` Black, `3` Hispanic, `4` Asian,
#'   `5` American Indian or Alaska Native.
#' @param age_menarche Age at menarche in years, or `"unknown"`.
#' @param age_first_birth Age at first live birth in years, `"nulliparous"`
#'   (no live births), or `"unknown"`.
#' @param n_biopsies Number of prior breast biopsies (`0`, `1`, `2` or more),
#'   or `"unknown"` (treated as `0`).
#' @param n_relatives Number of first-degree relatives (mother, sisters,
#'   daughters) with breast cancer (`0`, `1`, `2` or more), or `"unknown"`
#'   (treated as `0`).
#' @param atypical_hyperplasia Atypical hyperplasia on a prior biopsy:
#'   `"yes"`, `"no"`, or `"unknown"` (default). Ignored when `n_biopsies` is
#'   `0`/`"unknown"`.
#' @param hispanic_nativity Only used when `race` is Hispanic: `"us_born"`
#'   (default) or `"foreign_born"`.
#' @param asian_subgroup Only used when `race` is Asian: `"chinese"`,
#'   `"japanese"`, `"filipino"`, `"hawaiian"`, `"other_pacific_islander"`, or
#'   `"other_asian"` (default).
#'
#' @return A numeric vector of 5-year invasive breast cancer probabilities in
#'   `[0, 1]`, one element per (recycled) input row.
#'
#' @references
#' Gail MH, Brinton LA, Byar DP, et al. Projecting individualized
#' probabilities of developing breast cancer for white females who are being
#' examined annually. *J Natl Cancer Inst.* 1989;81(24):1879-1886.
#' \doi{10.1093/jnci/81.24.1879}
#'
#' Costantino JP, Gail MH, Pee D, et al. Validation studies for models
#' projecting the risk of invasive and total breast cancer incidence.
#' *J Natl Cancer Inst.* 1999;91(18):1541-1548. \doi{10.1093/jnci/91.18.1541}
#'
#' Gail MH, Costantino JP, Pee D, et al. Projecting individualized absolute
#' invasive breast cancer risk in African American women.
#' *J Natl Cancer Inst.* 2007;99(23):1782-1792. \doi{10.1093/jnci/djm223}
#'
#' Matsuno RK, Costantino JP, Ziegler RG, et al. Projecting individualized
#' absolute invasive breast cancer risk in Asian and Pacific Islander American
#' women. *J Natl Cancer Inst.* 2011;103(12):951-961. \doi{10.1093/jnci/djr154}
#'
#' Banegas MP, John EM, Slattery ML, et al. Projecting individualized absolute
#' invasive breast cancer risk in US Hispanic women. *J Natl Cancer Inst.*
#' 2017;109(2):djw215. \doi{10.1093/jnci/djw215}
#'
#' @examples
#' # A 50-year-old White woman, average risk factors
#' bcrat(age = 50, race = "White", age_menarche = 13, age_first_birth = 24,
#'       n_biopsies = 0, n_relatives = 0)
#'
#' # A 45-year-old foreign-born Hispanic woman
#' bcrat(age = 45, race = "Hispanic", hispanic_nativity = "foreign_born",
#'       age_menarche = 12, age_first_birth = 26, n_biopsies = 0,
#'       n_relatives = 1)
#' @export
bcrat <- function(age, race, age_menarche, age_first_birth, n_biopsies,
                  n_relatives, atypical_hyperplasia = "unknown",
                  hispanic_nativity = NULL, asian_subgroup = NULL) {

  n <- max(length(age), length(race), length(age_menarche), length(age_first_birth),
           length(n_biopsies), length(n_relatives), length(atypical_hyperplasia))
  rep_to_n <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) != n) stop("All inputs must have length 1 or a common length.", call. = FALSE)
    x
  }
  age                   <- as.numeric(rep_to_n(age))
  race                  <- rep_to_n(race)
  age_menarche          <- rep_to_n(age_menarche)
  age_first_birth       <- rep_to_n(age_first_birth)
  n_biopsies            <- rep_to_n(n_biopsies)
  n_relatives           <- rep_to_n(n_relatives)
  atypical_hyperplasia  <- rep_to_n(atypical_hyperplasia)
  hispanic_nativity     <- rep_to_n(hispanic_nativity)
  asian_subgroup        <- rep_to_n(asian_subgroup)

  # `atypical_hyperplasia` has a default, but a caller that passes it through
  # explicitly as NULL (e.g. a JSON payload with the field simply absent)
  # would otherwise bypass that default rather than trigger it.
  if (is.null(atypical_hyperplasia)) atypical_hyperplasia <- rep("unknown", n)

  if (any(age < 35 | age > 85 | age != round(age), na.rm = TRUE)) {
    stop("`age` must be a whole number from 35 to 85.", call. = FALSE)
  }

  group <- .bcrat_resolve_group(race, hispanic_nativity, asian_subgroup)

  nb <- .bcrat_recode_nb(n_biopsies, group)
  am <- .bcrat_recode_am(age_menarche, group)
  af <- .bcrat_recode_af(age_first_birth, group)
  nr <- .bcrat_recode_nr(n_relatives, group)
  r_hyp <- .bcrat_r_hyp(nb, atypical_hyperplasia)

  B <- .bcrat_beta[.bcrat_model_group(group), , drop = FALSE]
  lp1 <- nb * B[, "n_biopsies"] + am * B[, "age_menarche"] +
    af * B[, "age_first_birth"] + nr * B[, "n_relatives"] +
    af * nr * B[, "first_birth_x_relatives"] + log(r_hyp)
  lp2 <- lp1 + nb * B[, "age50_x_biopsies"]
  rr1 <- exp(lp1)
  rr2 <- exp(lp2)

  L1  <- .bcrat_lambda1[.bcrat_lambda_group(group), , drop = FALSE]
  L2  <- .bcrat_lambda2[.bcrat_lambda_group(group), , drop = FALSE]
  FAR <- .bcrat_far[.bcrat_model_group(group), , drop = FALSE]

  risk <- numeric(n)
  cum_hazard <- numeric(n)
  idx <- seq_len(n)
  for (k in 0:4) {
    year <- age + k
    band <- .bcrat_band(year)
    lambda1_y <- L1[cbind(idx, band)]
    lambda2_y <- L2[cbind(idx, band)]
    under50 <- year < 50
    far_y <- ifelse(under50, FAR[, 1], FAR[, 2])
    rr_y  <- ifelse(under50, rr1, rr2)
    hazard_y <- far_y * rr_y * lambda1_y
    total_y  <- hazard_y + lambda2_y
    risk <- risk + (hazard_y / total_y) * exp(-cum_hazard) * (1 - exp(-total_y))
    cum_hazard <- cum_hazard + total_y
  }
  unname(risk)
}
