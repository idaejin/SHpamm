shp_vcov_components <- function(object, method = c("V2", "V3")) {
  method <- match.arg(method)
  ped <- object$ped_data
  X <- stats::predict(object, newdata = ped, type = "lpmatrix")
  mu <- exp(drop(X %*% stats::coef(object)) + ped$offset)
  y <- ped$event
  w <- ped$weight
  n <- nrow(object$input_data)
  u <- shp_subject_rows(X * (w * (y - mu)), ped$id, n)
  u <- sweep(u, 2, colMeans(u))
  Vp <- object$Vp
  if (method == "V2") {
    return(list(covariance = Vp %*% crossprod(u) %*% Vp,
                psi = u, X = X, mu = mu, y = y))
  }
  if (object$weights_kind != "conditional") {
    stop("V3 requires a censoring PAMM fitted jointly by sh_pamm()",
         call. = FALSE)
  }
  if (object$weight_eval != "exact") {
    stop("V3 requires exact interval-averaged weights", call. = FALSE)
  }
  mc <- object$weights_model
  rebuilt <- shp_cond_weights(ped, mc, object$event_status,
                              eval_method = "exact", derivative = TRUE)
  if (max(abs(rebuilt$weight - w)) > 1e-7) {
    stop("Rebuilt weights do not agree with fitted weights", call. = FALSE)
  }
  cp <- object$censor_ped
  XC <- stats::predict(mc, newdata = cp, type = "lpmatrix")
  muC <- exp(drop(XC %*% stats::coef(mc)) + cp$offset)
  v <- shp_subject_rows(XC * (cp$event - muC), cp$id, n)
  v <- sweep(v, 2, colMeans(v))
  Ug <- crossprod(X * (y - mu), rebuilt$dw)
  psi <- u + v %*% mc$Vp %*% t(Ug)
  psi <- sweep(psi, 2, colMeans(psi))
  list(covariance = Vp %*% crossprod(psi) %*% Vp,
       psi = psi, v = v, dw = rebuilt$dw, X = X, mu = mu, y = y,
       U_gamma = Ug)
}

#' Variance of an SH-PAMM fit
#'
#' `V2` aggregates weighted Poisson scores by subject. `V3` also propagates
#' estimated censoring-PAMM coefficients through derivatives of exact weights.
#' For penalised fits with selected smoothness, both are approximations and do
#' not guarantee coverage of smooth functions. A full-refit bootstrap is not
#' returned as a coefficient covariance because resampled smooth bases need not
#' be identical; bootstrap predictions should be compared on a common grid.
#'
#' @param object A fitted `SHpamm` object.
#' @param method `"auto"` (V3 for internally fitted conditional weights, V2
#'   otherwise), `"V2"`, or `"V3"`. `"V4"` is reserved for a future release.
#' @param ... Unused.
#' @return Variance-covariance matrix for the model coefficients.
#' @examples
#' set.seed(4)
#' d <- data.frame(time = rexp(80), status = sample(0:2, 80, TRUE),
#'                 x = rnorm(80))
#' fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 5)
#' V <- vcov(fit, method = "V2")
#' dim(V)
#' @export
vcov.SHpamm <- function(object, method = c("auto", "V2", "V3", "V4"), ...) {
  method <- match.arg(method)
  if (method == "V4") {
    stop("V4 coefficient covariance is not available in SHpamm 0.1.0; refit and compare functionals on a common grid",
         call. = FALSE)
  }
  if (method == "auto") method <- if (object$weights_kind == "conditional" &&
                                    object$weight_eval == "exact") "V3" else "V2"
  V <- shp_vcov_components(object, method)$covariance
  dimnames(V) <- list(names(stats::coef(object)), names(stats::coef(object)))
  V
}
