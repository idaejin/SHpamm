shp_a <- function(lambda, duration) {
  x <- lambda * duration
  ans <- numeric(length(x))
  small <- abs(x) < 1e-7
  ans[small] <- duration[small] * (1 - x[small] / 2 + x[small]^2 / 6)
  ans[!small] <- -expm1(-x[!small]) / lambda[!small]
  ans
}

#' Compute exact interval-averaged weights for a piecewise exponential censoring model
#'
#' Implements the closed-form average of the retained competing-event IPCW
#' contribution. Rows must include all intervals from a subject's origin through
#' the analysis horizon and be ordered by subject and interval.
#'
#' @param ped Expanded subdistribution risk-set data from `shp_make_ped()`.
#' @param cens_fit A fitted `cens_pamm` object.
#' @param event_status Event-of-interest code.
#' @param eval_method `"exact"`, `"midpoint"`, or `"endpoint"`.
#' @param derivative Also return the derivative with respect to censoring GAM coefficients.
#' @return A list with `weight` and, when requested, matrix `dw`.
#' @keywords internal
shp_cond_weights <- function(ped, cens_fit, event_status,
                             eval_method = "exact", derivative = FALSE) {
  if (derivative && eval_method != "exact") {
    stop("Weight derivatives are implemented only for exact averaging", call. = FALSE)
  }
  X <- stats::predict(cens_fit, newdata = ped, type = "lpmatrix")
  lambda <- exp(drop(X %*% stats::coef(cens_fit)))
  if (any(!is.finite(lambda))) stop("Non-finite censoring hazard", call. = FALSE)
  D <- ped$tend - ped$tstart
  end <- ped$tstart + exp(ped$offset)
  Ti <- ped$.shp_time
  Hend <- ave(lambda * D, ped$id, FUN = cumsum)
  H0 <- Hend - lambda * D
  in_T <- Ti > ped$tstart & Ti <= ped$tend
  HT <- ave(ifelse(in_T, H0 + lambda * (Ti - ped$tstart), 0),
            ped$id, FUN = sum)
  comp <- ped$.shp_status != 0L & ped$.shp_status != event_status
  after <- comp & Ti <= ped$tstart
  straddle <- comp & Ti > ped$tstart & Ti < end
  w <- rep(1, nrow(ped))
  if (eval_method == "exact") {
    ia <- which(after)
    if (length(ia)) {
      w[ia] <- exp(-(H0[ia] - HT[ia])) * shp_a(lambda[ia], end[ia] - ped$tstart[ia]) /
        (end[ia] - ped$tstart[ia])
    }
    ib <- which(straddle)
    if (length(ib)) {
      d <- end[ib] - Ti[ib]
      w[ib] <- ((Ti[ib] - ped$tstart[ib]) + shp_a(lambda[ib], d)) /
        (end[ib] - ped$tstart[ib])
    }
  } else {
    at <- if (eval_method == "midpoint") (ped$tstart + end) / 2 else end
    i <- which(comp & Ti < at)
    if (length(i)) {
      w[i] <- exp(-(H0[i] + lambda[i] * (at[i] - ped$tstart[i]) - HT[i]))
    }
  }
  if (any(!is.finite(w)) || any(w < -1e-10 | w > 1 + 1e-10)) {
    stop("Invalid censoring weights", call. = FALSE)
  }
  w <- pmin(pmax(w, 0), 1)
  if (!derivative) return(list(weight = w, dw = NULL))

  q <- ncol(X)
  LX <- X * (lambda * D)
  Gend <- vapply(seq_len(q), function(j) ave(LX[, j], ped$id, FUN = cumsum),
                 numeric(nrow(ped)))
  Gend <- matrix(Gend, nrow = nrow(ped), ncol = q)
  G0 <- Gend - LX
  GTr <- (G0 + X * (lambda * (Ti - ped$tstart))) * in_T
  GT <- shp_subject_rows(GTr, ped$id, max(ped$id))[ped$id, , drop = FALSE]
  dw <- matrix(0, nrow(ped), q)
  ib <- which(straddle)
  if (length(ib)) {
    x <- end[ib] - Ti[ib]
    term <- (x * exp(-lambda[ib] * x) - shp_a(lambda[ib], x)) /
      (end[ib] - ped$tstart[ib])
    dw[ib, ] <- X[ib, , drop = FALSE] * term
  }
  ia <- which(after)
  if (length(ia)) {
    duration <- end[ia] - ped$tstart[ia]
    term <- exp(-(H0[ia] - HT[ia])) *
      (duration * exp(-lambda[ia] * duration) - shp_a(lambda[ia], duration)) /
      duration
    dw[ia, ] <- -(G0[ia, , drop = FALSE] - GT[ia, , drop = FALSE]) * w[ia] +
      X[ia, , drop = FALSE] * term
  }
  colnames(dw) <- colnames(X)
  list(weight = w, dw = dw)
}

shp_km_functions <- function(km_fit) {
  knots <- c(0, km_fit$time)
  height <- c(1, km_fit$surv)
  area <- c(0, cumsum(diff(knots) * head(height, -1L)))
  idx <- function(t) pmax(1L, findInterval(t, knots))
  list(
    surv = function(t) height[idx(t)],
    integral = function(t) {
      k <- idx(t)
      area[k] + height[k] * (t - knots[k])
    }
  )
}

shp_km_weights <- function(ped, km_fit, event_status,
                           eval_method = "exact") {
  fun <- shp_km_functions(km_fit)
  comp <- ped$.shp_status != 0L & ped$.shp_status != event_status
  w <- rep(1, nrow(ped))
  i <- which(comp)
  if (!length(i)) return(w)
  a <- ped$tstart[i]
  b <- a + exp(ped$offset[i])
  Ti <- ped$.shp_time[i]
  GT <- fun$surv(Ti)
  if (any(GT <= 0)) stop("Marginal censoring survival is zero at a competing event",
                         call. = FALSE)
  if (eval_method == "exact") {
    before <- pmax(0, pmin(Ti, b) - a)
    start <- pmax(a, Ti)
    after <- ifelse(b > start, (fun$integral(b) - fun$integral(start)) / GT, 0)
    w[i] <- (before + after) / (b - a)
  } else {
    at <- if (eval_method == "midpoint") (a + b) / 2 else b
    w[i] <- ifelse(at <= Ti, 1, fun$surv(at) / GT)
  }
  if (any(!is.finite(w)) || any(w < -1e-10 | w > 1 + 1e-10)) {
    stop("Invalid marginal censoring weights", call. = FALSE)
  }
  pmin(pmax(w, 0), 1)
}
