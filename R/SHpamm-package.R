#' Subdistribution Hazard Piecewise Exponential Additive Models
#'
#' `SHpamm` fits flexible Fine--Gray-type subdistribution hazard models using
#' piecewise exponential data and exact interval-averaged censoring weights.
#' It provides fitted CIFs and V2/V3 variance approximations. The software
#' handles baseline covariates and a common administrative horizon; it does not
#' currently handle delayed entry, time-dependent external covariates, or
#' subject-specific administrative closing dates.
#' @importFrom stats ave
#' @importFrom utils head tail
#'
#' @keywords internal
"_PACKAGE"
