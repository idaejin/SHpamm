#' Fit a piecewise exponential censoring model
#'
#' Fits a Poisson GAM for censoring before a known analysis horizon. Observed
#' events of every cause leave the censoring risk set without being counted as
#' censoring events. Observations at the common administrative horizon contribute
#' exposure but are not counted as random-censoring events.
#'
#' @param formula Right-hand-side formula for censoring predictors, for example
#'   `~ age + s(marker, bs = "ps")`. An optional left-hand side is ignored.
#' @param data One row per subject.
#' @param time_var,status_var Character names of time and integer status columns;
#'   status 0 denotes censoring and positive integers denote event causes.
#' @param event_status Positive code of the event of interest, used to set the
#'   default horizon.
#' @param tau Common analysis horizon. Defaults to the latest event of interest.
#' @param cut_points Number of cause-1 event-quantile intervals, or explicit
#'   increasing boundaries from 0 to `tau`.
#' @param smooth_method `"REML"`, `"GCV.Cp"`, or `"ML"`.
#' @param ... Additional arguments passed to [mgcv::bam()].
#' @return A `cens_pamm` object inheriting from `gam`. Its `ped_data` field holds
#'   the censoring PED used for fitting.
#' @examples
#' set.seed(2)
#' d <- data.frame(time = rexp(80), status = sample(0:2, 80, TRUE),
#'                 x = rnorm(80))
#' fit <- cens_pamm(~ x, d, "time", "status", tau = 1, cut_points = 5)
#' @export
cens_pamm <- function(formula, data, time_var, status_var, event_status = 1L,
                      tau = NULL, cut_points = 30,
                      smooth_method = c("REML", "GCV.Cp", "ML"), ...) {
  smooth_method <- match.arg(smooth_method)
  validated <- shp_validate_data(data, time_var, status_var, event_status, tau)
  data <- validated$data
  tau <- validated$tau
  cuts <- shp_cuts(cut_points, data, tau, event_status)
  if (!any(data$.shp_status == 0L & data$.shp_time < tau)) {
    stop("No pre-horizon censoring events: a conditional censoring PAMM cannot be fitted",
         call. = FALSE)
  }
  rhs <- shp_rhs(formula)
  covars <- intersect(all.vars(rhs), names(data))
  if (anyNA(data[covars])) stop("Missing censoring predictors are not supported",
                                  call. = FALSE)
  ped <- shp_make_ped(data, cuts, event_status, "censor")
  f <- shp_formula(rhs, "event", cuts)
  method <- if (smooth_method == "REML") "fREML" else smooth_method
  m <- mgcv::bam(f, data = ped, family = stats::poisson(),
                 method = method,
                 discrete = smooth_method == "REML" && shp_has_smooth(rhs, cuts), ...)
  if (!isTRUE(m$converged)) stop("Censoring PAMM did not converge", call. = FALSE)
  m$ped_data <- ped
  m$cut_points <- cuts
  m$tau <- tau
  m$time_var <- time_var
  m$status_var <- status_var
  m$event_status <- event_status
  class(m) <- c("cens_pamm", class(m))
  m
}
