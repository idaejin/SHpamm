#' Predict cumulative incidence from an SH-PAMM
#'
#' Integrates the fitted subdistribution hazard across the stored partition.
#' Pointwise intervals use a log-cumulative-hazard transformation and a
#' subject-level or stacked sandwich covariance. They are approximations for
#' penalised smooths, not simultaneous bands.
#'
#' @param object A fitted `SHpamm` object.
#' @param newdata Data frame of covariate profiles.
#' @param times Numeric times between 0 and the fitted horizon.
#' @param se.fit Include pointwise standard errors and intervals.
#' @param vcov_method `"auto"`, `"V2"`, or `"V3"`.
#' @param conf.level Confidence level between 0 and 1.
#' @param ... Unused.
#' @return Data frame containing the profiles, `time`, `cif`, and optionally
#'   `se`, `lower`, and `upper`.
#' @examples
#' set.seed(5)
#' d <- data.frame(time = rexp(90), status = sample(0:2, 90, TRUE),
#'                 x = rnorm(90))
#' fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 5)
#' predict_cif(fit, data.frame(x = 0), times = c(0.25, 0.5), se.fit = FALSE)
#' @export
predict_cif <- function(object, newdata, times, se.fit = TRUE,
                        vcov_method = c("auto", "V3", "V2"),
                        conf.level = 0.95, ...) {
  if (!inherits(object, "SHpamm")) stop("object must be an SHpamm fit", call. = FALSE)
  vcov_method <- match.arg(vcov_method)
  newdata <- as.data.frame(newdata)
  if (!nrow(newdata)) stop("newdata must contain at least one profile", call. = FALSE)
  if (!is.numeric(times) || !length(times) || anyNA(times) ||
      any(!is.finite(times) | times < 0 | times > object$tau)) {
    stop("times must be finite values between 0 and tau", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie between 0 and 1", call. = FALSE)
  }
  cuts <- object$cut_points
  J <- length(cuts) - 1L
  V <- if (se.fit) vcov.SHpamm(object, vcov_method) else NULL
  z <- if (se.fit) stats::qnorm((1 + conf.level) / 2) else NA_real_
  out <- vector("list", nrow(newdata) * length(times))
  index <- 0L
  for (i in seq_len(nrow(newdata))) {
    nd <- newdata[rep(i, J), , drop = FALSE]
    nd$tstart <- head(cuts, -1L)
    nd$tend <- tail(cuts, -1L)
    nd$interval <- factor(seq_len(J), levels = seq_len(J))
    nd$offset <- 0
    X <- stats::predict(object, newdata = nd, type = "lpmatrix")
    lambda <- exp(drop(X %*% stats::coef(object)))
    for (tt in times) {
      index <- index + 1L
      duration <- pmax(0, pmin(tt, tail(cuts, -1L)) - head(cuts, -1L))
      increments <- lambda * duration
      Lambda <- sum(increments)
      cif <- -expm1(-Lambda)
      row <- cbind(newdata[i, , drop = FALSE], time = tt, cif = cif)
      if (se.fit) {
        grad <- colSums(X * increments)
        seLambda <- sqrt(max(0, drop(crossprod(grad, V %*% grad))))
        se <- exp(-Lambda) * seLambda
        if (Lambda == 0) {
          lower <- upper <- 0
        } else {
          selog <- seLambda / Lambda
          lower <- -expm1(-Lambda * exp(-z * selog))
          upper <- -expm1(-Lambda * exp(z * selog))
        }
        row$se <- se
        row$lower <- lower
        row$upper <- upper
      }
      out[[index]] <- row
    }
  }
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}
