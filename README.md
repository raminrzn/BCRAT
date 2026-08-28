# BCRAT <img src="https://img.shields.io/badge/lifecycle-stable-brightgreen.svg" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/raminrzn/BCRAT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/raminrzn/BCRAT/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

A self-contained R implementation of the **NCI Breast Cancer Risk Assessment
Tool (BCRAT)**, commonly known as the **"Gail model,"** which predicts the
**5-year probability of invasive breast cancer**.

The package has **no dependencies** beyond base R and `jsonlite` (for the
ModelsCloud gateway), and exposes two layers:

1. **`bcrat()`** — a vectorised model function for interactive / batch use in R.
2. The **ModelsCloud API** (`model_run()`, `get_sample_input()`,
   `get_default_input()`) so the same package can be hosted as a prediction
   model on the [ModelsCloud](https://modelscloud.resp.core.ubc.ca/) platform
   (RESP Lab, UBC).

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("raminrzn/BCRAT")
```

---

## Quick start

```r
library(BCRAT)

# A 50-year-old White woman, average risk factors
bcrat(age = 50, race = "White", age_menarche = 13, age_first_birth = 24,
      n_biopsies = 0, n_relatives = 0)
#> [1] 0.0095...

# Score a whole cohort through the ModelsCloud API
model_run(get_sample_input())
#>   age     race age_menarche age_first_birth n_biopsies ... risk risk_percent
#> 1  50    White           13              24          0 ... 0.0x        x.xx
#> ...
```

---

## The model

BCRAT combines a **race-specific relative-risk model** with **race-specific
baseline hazards** (SEER invasive-breast-cancer incidence and NCHS competing,
non-breast-cancer mortality) in a single-year discrete competing-risk life
table (Gail et al., *JNCI* 1989, equation 1). For each of the five single
years of age from `age` to `age + 4`:

```
h   = F(t) * RR * lambda1        # attributable-risk-adjusted incidence hazard
tau = h + lambda2                # total hazard (breast cancer + competing)
pi  = (h / tau) * exp(-cumulative_tau) * (1 - exp(-tau))

5-year risk = sum(pi) over the five years
```

`lambda1`/`lambda2` are the race-specific baseline rates for the five-year age
band containing that year; `RR` is `RR1` (ages <50) or `RR2` (ages >=50); and
`F(t)` is one minus the attributable-risk fraction. The relative risk itself
is a small logistic model:

```
RR1 = exp(LP1)
RR2 = exp(LP1 + n_biopsies * age50_x_biopsies_beta)

LP1 = biopsies_beta   * n_biopsies_cat
    + menarche_beta   * age_menarche_cat
    + first_birth_beta * age_first_birth_cat
    + relatives_beta  * n_relatives_cat
    + first_birth_x_relatives_beta * age_first_birth_cat * n_relatives_cat
    + ln(hyperplasia_multiplier)
```

Continuous ages are first mapped to small categories (see **Input coding**
below); `hyperplasia_multiplier` is `0.93` (no), `1.82` (yes), or `1.00`
(unknown, or not applicable when there are no biopsies).

### Coefficients (relative-risk sub-model)

| Covariate | White | Black | Hispanic (US-born) | Hispanic (foreign-born) | Asian (all subgroups) |
|---|---:|---:|---:|---:|---:|
| n_biopsies | 0.5292641686 | 0.1822121131 | 0.0970783641 | 0.4798624017 | 0.55263612260619 |
| age_menarche | 0.0940103059 | 0.2672530336 | 0 | 0.2593922322 | 0.07499257592975 |
| age_first_birth | 0.2186262218 | 0 | 0.2318368334 | 0.4669246218 | 0.27638268294593 |
| n_relatives | 0.9583027845 | 0.4757242578 | 0.1666854410 | 0.9076679727 | 0.79185633720481 |
| age50 x n_biopsies | −0.2880424830 | −0.1119411682 | 0 | 0 | 0 |
| age_first_birth x n_relatives | −0.1908113865 | 0 | 0 | 0 | 0 |

**American Indian or Alaska Native** has no independently validated model —
NCI's own deployed tool substitutes the White coefficients (and White
baseline rates, below) in full, and this package does the same for exact
parity. **Black**'s coefficients being 0 for `age_first_birth` (and its
interaction) means that field never affects a Black estimate — the CARE study
model excludes it entirely. **US-born Hispanic**'s `age_menarche` coefficient
being 0 means that field never affects a US-born Hispanic estimate.

### Baseline rates (14 five-year age bands, [20,25) ... [85,90))

| Group | Incidence (`lambda1`) source | Mortality (`lambda2`) source | F(t), age&lt;50 / age&gt;=50 |
|---|---|---|---:|
| White | SEER White 1983-87 | NCHS White 1985-87 | 0.5788413 / 0.5788413 |
| Black | SEER Black 1994-98 | NCHS Black 1996-2000 | 0.72949880 / 0.74397137 |
| Hispanic (US-born) | CA Cancer Registry/SEER 1995-2004 | CA Cancer Registry/SEER 1995-2004 | 0.749294788397 / 0.778215491668 |
| Hispanic (foreign-born) | CA Cancer Registry/SEER 1995-2004 | CA Cancer Registry/SEER 1995-2004 | 0.428864989813 / 0.450352338746 |
| Asian (all 6 subgroups) | SEER-18 1998-2002 (per subgroup) | NCHS 1998-2002 (per subgroup) | 0.47519806426735 / 0.50316401683903 |
| American Indian/Alaska Native | = White | = White | = White |

The full 14-band numeric tables (and the six individual Asian subgroup
tables: Chinese, Japanese, Filipino, Hawaiian, Other Pacific Islander, Other
Asian) are in [`R/bcrat.R`](R/bcrat.R).

> **On accuracy and sourcing.** These coefficients and rate tables were
> cross-verified by independently executing NCI's own reference algorithm —
> distributed under GPL (>= 2) as the CRAN package
> ["BCRA"](https://cran.r-project.org/package=BCRA), maintained by the same
> person listed as BCRAT's technical-support contact on NCI's own site —
> against the published worked examples in the papers below and against
> [epiverse/bcra-js](https://github.com/epiverse/bcra-js)'s independent
> R-generated cross-validation fixtures. This package's R code is written
> from scratch against the published algorithm; it does not reuse BCRA's
> source.

---

## Race, ethnicity, and the categories that need one more field

`race` accepts `"White"`, `"Black"`, `"Hispanic"`, `"Asian"`, or
`"American Indian or Alaska Native"` (case-insensitive, with common aliases;
numeric codes `0`-`5` also work: `0`/`1` White, `2` Black, `3` Hispanic, `4`
Asian, `5` American Indian or Alaska Native). Two categories need one more
piece of information to reproduce NCI's own deployed tool exactly, because
the underlying science doesn't support a single pooled model for them —
**this package still works without it (sensible defaults apply), but the
result is only exactly right when the extra field is known and supplied:**

* **Hispanic** — US-born and foreign-born women have two entirely distinct
  fitted models (not one model with a nativity adjustment), confirmed
  directly from NCI's live production tool: selecting "Hispanic" requires a
  mandatory nativity choice with no combined/unspecified option. Set
  `hispanic_nativity = "us_born"` (the default here) or `"foreign_born"`.
* **Asian** — all six SEER subgroups share one relative-risk model but each
  has its own baseline incidence/mortality table, and NCI's tool has no
  pooled "combined Asian" option (confirmed the same way: a subgroup pick is
  mandatory in the live tool). Set `asian_subgroup` to `"chinese"`,
  `"japanese"`, `"filipino"`, `"hawaiian"`, `"other_pacific_islander"`, or
  `"other_asian"` (the default here — NCI's own catch-all for Asian
  ethnicities not individually tabulated).

`"Unknown"` race is mapped the same way NCI's own tool maps it: to the
American Indian/Alaska Native ("Other") bucket, i.e. the White model.

---

## Input coding

| Variable | Meaning | Coding |
|---|---|---|
| `age` | Age in years | whole number, 35-85 (the range NCI's own tool supports) |
| `race` | Race / ethnicity | label or numeric code (see above) |
| `age_menarche` | Age at menarche | years, or `"unknown"` |
| `age_first_birth` | Age at first live birth | years, `"nulliparous"`, or `"unknown"` |
| `n_biopsies` | Number of prior breast biopsies | `0`, `1`, `2` or more, or `"unknown"` |
| `n_relatives` | First-degree relatives (mother/sisters/daughters) with breast cancer | `0`, `1`, `2` or more, or `"unknown"` |
| `atypical_hyperplasia` | Atypical hyperplasia on a prior biopsy | `"yes"` · `"no"` · `"unknown"` (default); ignored if `n_biopsies` is 0 |
| `hispanic_nativity` | Only used when `race` is Hispanic | `"us_born"` (default) · `"foreign_born"` |
| `asian_subgroup` | Only used when `race` is Asian | `"chinese"` · `"japanese"` · `"filipino"` · `"hawaiian"` · `"other_pacific_islander"` · `"other_asian"` (default) |

A few coding nuances carried over faithfully from the official tool:
nulliparity is grouped with ages 25-29 (not the oldest age-at-first-birth
category); "unknown" values are always scored at the reference/lowest-risk
category; and Hispanic/Asian have their own category-pooling rules (e.g.
Hispanic pools "2 or more relatives" down to "1 or more").

---

## ModelsCloud entry points

The functions the ModelsCloud (pexa) executor calls — the package's hosted API:

| Function | Description |
|---|---|
| `model_run(model_input)` | Score a named list (one woman) or data frame (one row per woman); returns the input plus `risk` and `risk_percent`. |
| `get_sample_input(n)` | An example cohort, ready to pass to `model_run()`. |
| `get_default_input()` | One baseline woman to modify. |

Once hosted, end users call it through the
[`modelscloud`](https://github.com/resplab/modelscloud) client:

```r
library(modelscloud)
connect_to_model("raminrzn/bcrat", access_key = "YOUR_API_KEY")
result <- model_run(get_sample_input())
```

**Flexible input.** `model_run()` is forgiving about payload shape, so it
works with the `modelscloud` client *and* with raw API calls:

* Fields may be **wrapped** under `model_input` or passed as **top-level**
  named arguments (`do.call(model_run, funcInput)` style).
* NCI's own raw field names are understood as aliases: `T1`→`age`,
  `AgeMen`→`age_menarche`, `Age1st`→`age_first_birth`, `N_Biop`→`n_biopsies`,
  `N_Rels`→`n_relatives`, `HypPlas`→`atypical_hyperplasia`,
  `nativity`→`hispanic_nativity`, `ethnicity`/`subgroup`→`asian_subgroup`.
* Unknown extra fields are ignored.

### Raw HTTP

The ModelsCloud executor runs `do.call(model_run, funcInput)`, so wrap the
fields under `model_input`:

```bash
curl -X POST https://core.modelscloud.resp.core.ubc.ca/call/v2/<ns>/bcrat \
  -H "Authorization: Bearer <ACCESS_KEY_TOKEN>" \
  -H "Content-Type: application/json" \
  -d '{"funcInput": {"model_input": {
        "age": 50, "race": "hispanic", "age_menarche": 12,
        "age_first_birth": 26, "n_biopsies": 0, "n_relatives": 1,
        "hispanic_nativity": "foreign_born"
      }}}'
```

---

## Clinical interpretation

BCRAT estimates 5-year invasive breast cancer risk; a commonly used
chemoprevention-eligibility threshold in US clinical guidance is **>=1.67%**.
This package returns the raw probability and does **not** impose a threshold
— interpretation is left to the user / program. NCI's own documentation
notes that Hispanic estimates in particular "are subject to greater
uncertainty than those for White and Black/African American women," and that
American Indian/Alaska Native estimates are the White model's, substituted
for lack of an independently validated alternative.

> This software is for research use. It is **not** a medical device and is not a
> substitute for clinical judgement.

---

## References

> Gail MH, Brinton LA, Byar DP, et al. Projecting individualized
> probabilities of developing breast cancer for white females who are being
> examined annually. *J Natl Cancer Inst.* 1989;81(24):1879–1886.
> doi:[10.1093/jnci/81.24.1879](https://doi.org/10.1093/jnci/81.24.1879)

> Costantino JP, Gail MH, Pee D, et al. Validation studies for models
> projecting the risk of invasive and total breast cancer incidence.
> *J Natl Cancer Inst.* 1999;91(18):1541–1548.
> doi:[10.1093/jnci/91.18.1541](https://doi.org/10.1093/jnci/91.18.1541)

> Gail MH, Costantino JP, Pee D, et al. Projecting individualized absolute
> invasive breast cancer risk in African American women.
> *J Natl Cancer Inst.* 2007;99(23):1782–1792.
> doi:[10.1093/jnci/djm223](https://doi.org/10.1093/jnci/djm223)

> Matsuno RK, Costantino JP, Ziegler RG, et al. Projecting individualized
> absolute invasive breast cancer risk in Asian and Pacific Islander American
> women. *J Natl Cancer Inst.* 2011;103(12):951–961.
> doi:[10.1093/jnci/djr154](https://doi.org/10.1093/jnci/djr154)

> Banegas MP, John EM, Slattery ML, et al. Projecting individualized absolute
> invasive breast cancer risk in US Hispanic women. *J Natl Cancer Inst.*
> 2017;109(2):djw215.
> doi:[10.1093/jnci/djw215](https://doi.org/10.1093/jnci/djw215)

## License

GPL-3. Model © its original authors; package implementation © Ramin Rezaeianzadeh.
