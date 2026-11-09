toy_competing <- function(n = 120L, seed = 9L) {
  set.seed(seed)
  x <- rnorm(n)
  t1 <- rexp(n, .4 * exp(.4 * x))
  t2 <- rexp(n, .28)
  cens <- rexp(n, .22 * exp(.2 * x))
  time <- pmin(t1, t2, cens, 1.5)
  status <- ifelse(time == t1, 1L, ifelse(time == t2, 2L, 0L))
  data.frame(time, status, x)
}
