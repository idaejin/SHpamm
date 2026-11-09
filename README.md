# SHpamm

`SHpamm` fits flexible subdistribution hazard models for competing risks using a weighted piecewise exponential additive model. It provides exact interval-averaged censoring weights, `mgcv` smooths, cumulative incidence predictions, and V2/V3 sandwich variance approximations.

This is a development release. It has not been submitted to CRAN. The associated manuscript is under preparation; no journal acceptance is claimed.

## Install

```r
install.packages("remotes")
remotes::install_github("idaejin/SHpamm")
```

## Data format and first fit

Supply one row per subject. The time column must be positive. In the status column, `0` means censored, `1` is the event of interest, and other positive integers are competing events. **Keep censored subjects in the data.** Set a clinically justified common horizon `tau` and identify known administrative observations separately from earlier censoring in the source data.

```r
library(SHpamm)

set.seed(42)
n <- 160
x <- rnorm(n)
t1 <- rexp(n, 0.30 * exp(0.4 * x))
t2 <- rexp(n, 0.22)
cens <- rexp(n, 0.16)
time <- pmin(t1, t2, cens, 2)
status <- ifelse(time == t1, 1L, ifelse(time == t2, 2L, 0L))
d <- data.frame(time, status, x)

fit <- sh_pamm(
  ~ x,
  data = d, time_var = "time", status_var = "status",
  event_status = 1L, tau = 2, cut_points = 10,
  weights_from = "conditional", cens_formula = ~ x
)
predict_cif(fit, data.frame(x = c(-1, 0, 1)), times = c(1, 2))
summary(censor_support(fit))
```

The SH model accepts `mgcv` terms such as `s(x, bs="ps")`; a time interaction can use the expanded-data endpoint `tend`, for example `s(tend, by=treatment, bs="ps")`. The time baseline is added automatically. A supplied `Surv()` response is ignored: `time_var` and the original integer `status_var` determine the competing-risk data.

`weights_from = NULL` uses a marginal Kaplan–Meier censoring estimate. `weights_from = "conditional"` fits a censoring PAMM from `cens_formula`; an existing `cens_pamm` fit may also be supplied for point estimation. V3 is available only when the conditional censoring model was fitted inside `sh_pamm()` with exact weights. V2 treats the fitted weights as known. Both are approximations for penalised smooths. `vcov(..., method="V4")` is intentionally unavailable in version 0.1.0 because resampled smooth bases do not generally share a coefficient parametrisation.

The package assumes conditionally independent censoring, adequate support, baseline covariates, and a common known analysis horizon. It does not handle delayed entry, external time-dependent covariates, or unknown subject-specific administrative closing dates. `censor_support()` reports fitted survival only; it cannot prove population positivity. The package provides an approximate PSH score diagnostic through `sh_diagnostic()` for exact conditional-weight fits.

## Rotterdam data recipe

The public `survival::rotterdam` data have separate recurrence and death times. A usable competing-risk status must be constructed before fitting; `Surv(dtime, status == 1)` and removal of censored rows would be incorrect. Package examples use simulated data so that the status definition is explicit.

For the five-year endpoint in the accompanying manuscript, the analysis data can be prepared as follows. Here `late` identifies records with death reported after the recorded recurrence follow-up time and no reported recurrence; those records leave follow-up at `rtime`.

```r
r <- survival::rotterdam
late <- r$recur == 0 & r$death == 1 & r$dtime > r$rtime
rot <- data.frame(
  time = ifelse(r$recur == 1 | late, r$rtime, r$dtime) / 365.25,
  status = ifelse(r$recur == 1, 1L, 2L * r$death * !late),
  age = r$age,
  nodes = pmin(r$nodes, 20)
)
table(rot$status)
rot_fit <- sh_pamm(~ age + nodes, rot, "time", "status",
                   tau = 5, cut_points = 20,
                   weights_from = "conditional",
                   cens_formula = ~ age + nodes)
predict_cif(rot_fit, data.frame(age = 55, nodes = 5), times = 5)
```

## Development

Run `R CMD build SHpamm` and `R CMD check SHpamm_0.1.0.tar.gz` from the parent directory. See `vignette("intro", package="SHpamm")` for a worked example.
