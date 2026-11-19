for (f in list.files(file.path("..", "..", "R"), full.names = TRUE)) source(f)

test_that("log_returns refuses non-positive prices unless told otherwise", {
  expect_error(log_returns(c(10, -1, 12)), "non-positive")
  r <- suppressMessages(log_returns(c(10, -1, 12, 12), on_nonpositive = "na"))
  expect_true(all(is.na(r[1:3])))
  expect_equal(r[4], 0)
  expect_equal(log_returns(c(100, 110))[2], 100 * log(1.1))
})

test_that("Student-t log-likelihood matches dt()", {
  set.seed(1); r <- rt(500, 5) * 2 + 0.3
  par <- c(0.3, log(2), log(5 - 2))
  expect_equal(-t_nll(par, r), sum(dt((r - 0.3) / 2, 5, log = TRUE) - log(2)))
})

test_that("unit-variance t density integrates to 1 and has variance 1", {
  f <- function(z) exp(std_t_logdens(z, 6))
  expect_equal(integrate(f, -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(function(z) z^2 * f(z), -Inf, Inf)$value, 1, tolerance = 1e-6)
})

test_that("t expected shortfall formula agrees with Monte Carlo", {
  set.seed(2); x <- rt(2e6, 4); q <- qt(0.01, 4)
  es_formula <- -dt(q, 4) / 0.01 * (4 + q^2) / 3
  expect_equal(es_formula, mean(x[x <= q]), tolerance = 0.02)
})

test_that("GARCH(1,1)-t MLE recovers simulated parameters", {
  truth <- c(omega = 0.05, alpha = 0.08, beta = 0.90, nu = 6)
  r <- simulate_garch_t(8000, truth[1], truth[2], truth[3], truth[4], seed = 3)
  g <- fit_garch_t(r)
  expect_equal(unname(g$coef["alpha"]), 0.08, tolerance = 0.25)
  expect_equal(unname(g$coef["beta"]), 0.90, tolerance = 0.03)
  expect_equal(unname(g$coef["nu"]), 6, tolerance = 0.3)
  expect_true(all(is.finite(g$se)))
})

test_that("Kupiec test does not reject exact coverage and handles zero hits", {
  expect_gt(kupiec_test(10, 1000, 0.01)$p_value, 0.99)
  expect_true(is.finite(kupiec_test(0, 1000, 0.01)$LR))
  expect_lt(kupiec_test(30, 1000, 0.01)$p_value, 0.001)
})

test_that("weekend events map to the next trading day", {
  td <- as.Date(c("2019-09-12", "2019-09-13", "2019-09-16", "2019-09-17"))
  expect_equal(map_day0(as.Date(c("2019-09-14", "2019-09-13", "2019-09-18")), td), c(3, 2, NA))
})

test_that("overlapping events are flagged", {
  expect_equal(flag_overlap(c(100, 103, 200, 300, 301)), c(FALSE, TRUE, FALSE, FALSE, TRUE))
})

test_that("AR(1) half-life is recovered on simulated data", {
  set.seed(4); x <- as.numeric(arima.sim(list(ar = 0.95), n = 20000))
  res <- ar1_half_life(x)
  expect_equal(res$half_life, log(0.5) / log(0.95), tolerance = 0.15)
})
