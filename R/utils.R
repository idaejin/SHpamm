shp_rhs <- function(formula) {
  if (!inherits(formula, "formula")) stop("formula must be a formula", call. = FALSE)
  rhs <- if (length(formula) == 3L) formula[[3L]] else formula[[2L]]
  stats::as.formula(call("~", rhs), env = environment(formula))
}

shp_formula <- function(rhs, response, cuts) {
  k <- min(8L, length(cuts) - 1L)
  baseline <- if (k >= 4L) sprintf("s(tend, bs='ps', k=%d)", k) else "interval"
  rhs_text <- paste(deparse(rhs[[2L]], width.cutoff = 500L), collapse = "")
  f <- stats::as.formula(paste(response, "~", baseline, "+", rhs_text,
                               "+ offset(offset)"), env = environment(rhs))
  environment(f) <- list2env(list(s = mgcv::s, te = mgcv::te,
                                  ti = mgcv::ti, t2 = mgcv::t2),
                             parent = environment(rhs))
  f
}

shp_has_smooth <- function(rhs, cuts) {
  length(cuts) - 1L >= 4L ||
    grepl("(s|te|ti|t2)\\(", paste(deparse(rhs), collapse = ""))
}

shp_validate_data <- function(data, time_var, status_var, event_status, tau) {
  data <- as.data.frame(data)
  if (!is.character(time_var) || length(time_var) != 1L ||
      !is.character(status_var) || length(status_var) != 1L ||
      !all(c(time_var, status_var) %in% names(data))) {
    stop("time_var and status_var must name columns in data", call. = FALSE)
  }
  reserved <- c("id", "tstart", "tend", "interval", "offset", "event",
                "weight", ".shp_time", ".shp_status")
  if (any(reserved %in% names(data))) {
    stop("Input data use reserved package column names: ",
         paste(intersect(reserved, names(data)), collapse = ", "), call. = FALSE)
  }
  time <- data[[time_var]]
  status <- data[[status_var]]
  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time) | time <= 0)) {
    stop("Event times must be finite and strictly positive", call. = FALSE)
  }
  if (!is.numeric(status) || anyNA(status) || any(status < 0 | status != floor(status))) {
    stop("status_var must contain integer codes: 0 censoring, positive causes",
         call. = FALSE)
  }
  if (!is.numeric(event_status) || length(event_status) != 1L ||
      !is.finite(event_status) || event_status <= 0 || event_status != floor(event_status)) {
    stop("event_status must be a positive integer cause code", call. = FALSE)
  }
  if (!any(status == event_status)) stop("No event of interest was observed", call. = FALSE)
  if (is.null(tau)) tau <- max(time[status == event_status])
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be a positive finite number", call. = FALSE)
  }
  data$.shp_time <- pmin(time, tau)
  data$.shp_status <- ifelse(time > tau, 0L, as.integer(status))
  if (!any(data$.shp_status == event_status)) {
    stop("No event of interest occurs within tau", call. = FALSE)
  }
  list(data = data, tau = tau)
}

shp_cuts <- function(cut_points, data, tau, event_status) {
  if (!is.numeric(cut_points) || anyNA(cut_points) ||
      any(!is.finite(cut_points))) stop("cut_points must be finite numeric values", call. = FALSE)
  if (length(cut_points) == 1L) {
    j <- as.integer(cut_points)
    if (j < 1L || j != cut_points) stop("cut_points count must be a positive integer",
                                      call. = FALSE)
    ev <- data$.shp_time[data$.shp_status == event_status & data$.shp_time < tau]
    interior <- if (length(ev)) stats::quantile(ev, probs = seq(0, 1, length.out = j + 1L)[-c(1L, j + 1L)],
                                                names = FALSE, type = 8) else numeric()
    cuts <- sort(unique(c(0, interior[interior > 0 & interior < tau], tau)))
  } else {
    cuts <- as.numeric(cut_points)
    if (length(cuts) < 2L || is.unsorted(cuts, strictly = TRUE) ||
        cuts[1L] != 0 || abs(tail(cuts, 1L) - tau) > 1e-10) {
      stop("Explicit cut_points must increase strictly from 0 to tau", call. = FALSE)
    }
  }
  if (length(cuts) < 2L) stop("Partition has no interval", call. = FALSE)
  cuts
}

shp_subject_rows <- function(x, id, n) {
  ans <- rowsum(x, id, reorder = TRUE)
  out <- matrix(0, n, ncol(ans), dimnames = list(as.character(seq_len(n)), colnames(ans)))
  out[as.integer(rownames(ans)), ] <- ans
  out
}

shp_wald <- function(score, covariance, tol = 1e-8) {
  e <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  keep <- e$values > max(e$values) * tol & e$values > 0
  if (!any(keep)) stop("Score covariance has no positive rank", call. = FALSE)
  z <- drop(crossprod(e$vectors[, keep, drop = FALSE], score))
  stat <- sum(z^2 / e$values[keep])
  list(statistic = stat, parameter = sum(keep),
       p.value = stats::pchisq(stat, sum(keep), lower.tail = FALSE))
}
