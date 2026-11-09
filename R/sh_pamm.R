#' Fit a subdistribution hazard PAMM
#'
#' Fits a weighted Poisson GAM on an expanded subdistribution risk set. Subjects
#' with a competing event remain at risk to `tau`; their exposure is weighted by
#' the exact interval average of the conditional or marginal censoring-survival
#' ratio. The formula describes the hazard covariates on the right-hand side;
#' a smooth time baseline is added automatically. An optional response such as
#' `Surv(time, status == 1)` is ignored: the original integer `status_var` is
#' always required so that censoring and competing events remain distinct.
#'
#' @param formula Right-hand-side model formula, e.g. `~ age + s(marker, bs="ps")`.
#'   Time interactions may use the PED endpoint `tend`, e.g.
#'   `s(tend, by = treatment, bs = "ps")`.
#' @param data One row per subject, including censored subjects.
#' @param time_var,status_var Character names of time and integer status columns;
#'   0 is censoring and positive codes are event causes.
#' @param event_status Code of the event of interest (default 1).
#' @param tau Common analysis horizon; defaults to the latest event of interest.
#' @param cut_points Number of event-quantile intervals (default 30), or explicit
#'   boundaries starting at 0 and ending at `tau`.
#' @param weights_from `NULL` for marginal Kaplan–Meier weights,
#'   `"conditional"` to fit a censoring PAMM, or an existing `cens_pamm` object.
#' @param cens_formula Censoring predictors for `weights_from="conditional"`.
#' @param weight_eval `"exact"` (default), `"midpoint"`, or `"endpoint"`.
#' @param smooth_method Smoothing selection: `"REML"`, `"GCV.Cp"`, or `"ML"`.
#' @param ... Additional arguments passed to [mgcv::bam()] for the SH model.
#' @return An `SHpamm` object inheriting from `gam`, with `ped_data`,
#'   `weights_model`, `cut_points`, and `tau` fields.
#' @examples
#' set.seed(7)
#' d <- data.frame(time = rexp(100), status = sample(0:2, 100, TRUE),
#'                 x = rnorm(100))
#' fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 6)
#' predict_cif(fit, data.frame(x = 0), times = c(0.25, 0.5), se.fit = FALSE)
#' @export
sh_pamm <- function(formula, data, time_var, status_var, event_status = 1L,
                    tau = NULL, cut_points = 30, weights_from = NULL,
                    cens_formula = NULL,
                    weight_eval = c("exact", "midpoint", "endpoint"),
                    smooth_method = c("REML", "GCV.Cp", "ML"), ...) {
  call <- match.call()
  weight_eval <- match.arg(weight_eval)
  smooth_method <- match.arg(smooth_method)
  raw_data <- as.data.frame(data)
  validated <- shp_validate_data(data, time_var, status_var, event_status, tau)
  data <- validated$data
  tau <- validated$tau
  cuts <- shp_cuts(cut_points, data, tau, event_status)
  rhs <- shp_rhs(formula)
  covars <- intersect(all.vars(rhs), names(data))
  if (anyNA(data[covars])) stop("Missing SH predictors are not supported", call. = FALSE)
  ped <- shp_make_ped(data, cuts, event_status, "sh")

  if (is.null(weights_from)) {
    censor_event <- data$.shp_status == 0L & data$.shp_time < tau
    weights_model <- survival::survfit(survival::Surv(data$.shp_time,
                                                    censor_event) ~ 1)
    ped$weight <- shp_km_weights(ped, weights_model, event_status, weight_eval)
    weights_kind <- "marginal"
  } else if (identical(weights_from, "conditional")) {
    if (is.null(cens_formula)) stop("cens_formula is required for conditional weights",
                                     call. = FALSE)
    weights_model <- cens_pamm(cens_formula, raw_data, time_var, status_var,
                               event_status, tau, cuts, smooth_method)
    ped$weight <- shp_cond_weights(ped, weights_model, event_status,
                                   weight_eval)$weight
    weights_kind <- "conditional"
  } else if (inherits(weights_from, "cens_pamm")) {
    if (length(weights_from$cut_points) != length(cuts) ||
        any(abs(weights_from$cut_points - cuts) > 1e-10)) {
      stop("Pre-fitted censoring model must use the same cut points", call. = FALSE)
    }
    weights_model <- weights_from
    ped$weight <- shp_cond_weights(ped, weights_model, event_status,
                                   weight_eval)$weight
    weights_kind <- "prefitted"
  } else stop("weights_from must be NULL, 'conditional', or a cens_pamm object",
              call. = FALSE)

  f <- shp_formula(rhs, "event", cuts)
  method <- if (smooth_method == "REML") "fREML" else smooth_method
  weights_vec <- ped$weight
  environment(f)$weights_vec <- weights_vec
  fit <- mgcv::bam(f, data = ped, family = stats::poisson(),
                   weights = weights_vec, method = method,
                   discrete = smooth_method == "REML" && shp_has_smooth(rhs, cuts), ...)
  if (!isTRUE(fit$converged)) stop("SH-PAMM did not converge", call. = FALSE)
  fit$call <- call
  fit$ped_data <- ped
  fit$weights_model <- weights_model
  fit$weights_kind <- weights_kind
  fit$censor_ped <- if (inherits(weights_model, "cens_pamm"))
    weights_model$ped_data else NULL
  fit$cut_points <- cuts
  fit$tau <- tau
  fit$event_status <- event_status
  fit$weight_eval <- weight_eval
  fit$smooth_method <- smooth_method
  fit$input_data <- data
  class(fit) <- c("SHpamm", class(fit))
  fit
}
