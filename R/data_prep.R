#' Build the expanded subdistribution risk set
#'
#' Subjects with a competing event remain in the subdistribution risk set to
#' `tau`; subjects with the event of interest or earlier censoring leave at their
#' observed time. A known common administrative boundary is exposure, not an
#' event in the censoring model.
#'
#' @param data One row per subject, prepared by `shp_validate_data()`.
#' @param cuts Partition including 0 and the analysis horizon.
#' @param event_status Positive integer code of the event of interest.
#' @param mode Either `"sh"` or `"censor"`.
#' @return Expanded data frame with exact exposures and event indicators.
#' @keywords internal
shp_make_ped <- function(data, cuts, event_status, mode = c("sh", "censor")) {
  mode <- match.arg(mode)
  n <- nrow(data)
  j <- length(cuts) - 1L
  grid <- data.table::data.table(id = rep(seq_len(n), each = j),
                                 .j = rep(seq_len(j), times = n))
  i <- grid$id
  k <- grid$.j
  end_time <- if (mode == "sh") {
    ifelse(data$.shp_status != 0L & data$.shp_status != event_status,
           tail(cuts, 1L), data$.shp_time)
  } else data$.shp_time
  exposure <- pmin(end_time[i], cuts[k + 1L]) - cuts[k]
  keep <- exposure > 0
  i <- i[keep]
  k <- k[keep]
  exposure <- exposure[keep]
  ped <- data[i, , drop = FALSE]
  rownames(ped) <- NULL
  ped$id <- i
  ped$tstart <- cuts[k]
  ped$tend <- cuts[k + 1L]
  ped$interval <- factor(k, levels = seq_len(j))
  ped$offset <- log(exposure)
  is_event <- if (mode == "sh") data$.shp_status[i] == event_status else
    data$.shp_status[i] == 0L & data$.shp_time[i] < tail(cuts, 1L)
  ped$event <- as.integer(is_event & data$.shp_time[i] > ped$tstart &
                            data$.shp_time[i] <= ped$tend)
  ped$weight <- 1
  expected <- if (mode == "sh") sum(data$.shp_status == event_status) else
    sum(data$.shp_status == 0L & data$.shp_time < tail(cuts, 1L))
  if (sum(ped$event) != expected) stop("PED event count does not match subject data",
                                         call. = FALSE)
  ped
}
