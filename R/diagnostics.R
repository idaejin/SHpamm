#' Fitted censoring support at the analysis horizon
#'
#' Evaluates fitted conditional censoring survival at `tau` over the observed
#' covariate profiles. This is a model diagnostic and cannot establish
#' population positivity or conditional independent censoring.
#'
#' @param object A fitted `SHpamm` object with conditional censoring weights.
#' @return Numeric fitted survival probabilities, one per input subject.
#' @examples
#' set.seed(3)
#' d <- data.frame(time = rexp(100), status = sample(0:2, 100, TRUE),
#'                 x = rnorm(100))
#' fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 5,
#'                weights_from = "conditional", cens_formula = ~ x)
#' summary(censor_support(fit))
#' @export
censor_support <- function(object) {
  if (!inherits(object, "SHpamm") ||
      !inherits(object$weights_model, "cens_pamm")) {
    stop("A conditional censoring PAMM is required", call. = FALSE)
  }
  cuts <- object$cut_points
  J <- length(cuts) - 1L
  n <- nrow(object$input_data)
  nd <- object$input_data[rep(seq_len(n), each = J), , drop = FALSE]
  nd$tstart <- rep(head(cuts, -1L), n)
  nd$tend <- rep(tail(cuts, -1L), n)
  nd$interval <- factor(rep(seq_len(J), n), levels = seq_len(J))
  nd$offset <- 0
  h <- exp(drop(stats::predict(object$weights_model, newdata = nd,
                               type = "lpmatrix") %*% stats::coef(object$weights_model)))
  H <- rowsum(h * (nd$tend - nd$tstart), rep(seq_len(n), each = J),
              reorder = TRUE)
  exp(-drop(H))
}

#' Score diagnostic for a time-varying subdistribution effect
#'
#' Tests an added fixed-basis P-spline interaction with time after fitting a
#' null model without that interaction. The added interaction is not assigned a
#' smoothing parameter. The statistic uses the stacked censoring correction;
#' nuisance smoothness selected in the null fit makes its chi-squared reference
#' approximate. The test requires exact conditional weights and does not protect
#' against misspecified SH or censoring models.
#'
#' @param object Null `SHpamm` fit with an internally fitted conditional
#'   censoring PAMM and no time interaction for `var`.
#' @param var Name of a numeric baseline covariate to test.
#' @param df Requested interaction basis rank (default 5).
#' @return An `htest` object with chi-squared statistic, effective rank and
#'   approximate p-value.
#' @examples
#' set.seed(12)
#' d <- data.frame(time = rexp(100), status = sample(0:2, 100, TRUE),
#'                 x = rnorm(100))
#' fit <- sh_pamm(~ x, d, "time", "status", tau = 1, cut_points = 6,
#'                weights_from = "conditional", cens_formula = ~ x)
#' sh_diagnostic(fit, "x")
#' @export
sh_diagnostic <- function(object, var, df = 5L) {
  if (!inherits(object, "SHpamm") || object$weights_kind != "conditional" ||
      object$weight_eval != "exact") {
    stop("Score diagnostic requires exact, internally fitted conditional weights",
         call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L ||
      !var %in% names(object$input_data) ||
      !is.numeric(object$input_data[[var]])) {
    stop("var must name a numeric baseline covariate", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || df < 2L || df != floor(df)) {
    stop("df must be an integer of at least two", call. = FALSE)
  }
  labels <- vapply(object$smooth, function(z) z$label, character(1))
  if (any(grepl(paste0("by=", var), gsub(" ", "", labels), fixed = TRUE))) {
    stop("The null fit already contains a time interaction for var", call. = FALSE)
  }
  ped <- object$ped_data
  spec <- eval(substitute(mgcv::s(tend, by = V, bs = "ps", k = K, pc = 0),
                          list(V = as.name(var), K = as.integer(df + 1L))))
  B <- mgcv::smoothCon(spec, data = ped, absorb.cons = TRUE)[[1L]]$X
  comp <- shp_vcov_components(object, "V3")
  X <- comp$X
  y <- comp$y
  mu <- comp$mu
  w <- ped$weight
  n <- nrow(object$input_data)
  ua <- shp_subject_rows(B * (w * (y - mu)), ped$id, n)
  ua <- sweep(ua, 2, colMeans(ua))
  theta_score <- colSums(X * (w * (y - mu)))
  A <- crossprod(B, X * (w * mu)) %*% object$Vp
  score <- colSums(B * (w * (y - mu))) - drop(A %*% theta_score)
  Uag <- crossprod(B * (y - mu), comp$dw)
  psi_a <- ua + comp$v %*% object$weights_model$Vp %*% t(Uag)
  eff <- psi_a - comp$psi %*% t(A)
  eff <- sweep(eff, 2, colMeans(eff))
  test <- shp_wald(score, crossprod(eff))
  structure(list(statistic = c(`X-squared` = test$statistic),
                 parameter = c(df = test$parameter), p.value = test$p.value,
                 method = "Approximate fixed-basis SH time-interaction score diagnostic",
                 data.name = var), class = "htest")
}
