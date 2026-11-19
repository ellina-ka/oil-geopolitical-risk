# Returns, return distributions and tail risk.

#' Log returns in percent between consecutive observations.
#' Non-positive prices make log returns undefined (WTI settled at -37.63 on
#' 2020-04-20). Default is to stop; on_nonpositive = "na" sets every return
#' that touches a non-positive price to NA and reports how many were dropped.
log_returns <- function(price, scale = 100, on_nonpositive = c("error", "na")) {
  on_nonpositive <- match.arg(on_nonpositive)
  if (!is.numeric(price)) stop("price must be numeric")
  bad <- !is.na(price) & price <= 0
  if (any(bad)) {
    if (on_nonpositive == "error") {
      stop(sum(bad), " non-positive price(s): log returns undefined. ",
           "Use on_nonpositive = 'na' to drop them explicitly.")
    }
    message(sum(bad), " non-positive price(s) set to NA; adjacent returns dropped.")
    price[bad] <- NA
  }
  c(NA_real_, scale * diff(log(price)))
}

#' Sample moments and Jarque-Bera normality test (population moments).
describe_returns <- function(r) {
  r <- r[is.finite(r)]
  n <- length(r)
  if (n < 30) stop("need at least 30 finite returns, got ", n)
  m <- mean(r)
  s <- sqrt(mean((r - m)^2))
  skew <- mean((r - m)^3) / s^3
  exkurt <- mean((r - m)^4) / s^4 - 3
  jb <- n / 6 * (skew^2 + exkurt^2 / 4)
  data.frame(n = n, mean = m, sd = sd(r), skewness = skew,
             excess_kurtosis = exkurt, jarque_bera = jb,
             jb_p_value = stats::pchisq(jb, df = 2, lower.tail = FALSE),
             min = min(r), max = max(r))
}

fit_normal <- function(r) {
  r <- r[is.finite(r)]
  mu <- mean(r)
  s <- sqrt(mean((r - mu)^2)) # MLE
  ll <- sum(stats::dnorm(r, mu, s, log = TRUE))
  list(mu = mu, sigma = s, loglik = ll, aic = 2 * 2 - 2 * ll, n = length(r))
}

# Location-scale Student-t. nu = 2 + exp(.) keeps the variance finite.
t_nll <- function(par, r) {
  mu <- par[1]; s <- exp(par[2]); nu <- 2 + exp(par[3])
  -sum(stats::dt((r - mu) / s, df = nu, log = TRUE) - log(s))
}

fit_t <- function(r) {
  r <- r[is.finite(r)]
  start <- c(stats::median(r), log(stats::mad(r)), log(3))
  opt <- stats::optim(start, t_nll, r = r, method = "BFGS",
                      control = list(maxit = 500))
  if (opt$convergence != 0) warning("fit_t: optim did not converge (code ", opt$convergence, ")")
  list(mu = opt$par[1], scale = exp(opt$par[2]), nu = 2 + exp(opt$par[3]),
       loglik = -opt$value, aic = 2 * 3 + 2 * opt$value, n = length(r))
}

#' One-day Value-at-Risk and Expected Shortfall at tail probability `p`,
#' reported as positive losses (in return units).
var_es <- function(r, p = 0.01) {
  stopifnot(p > 0, p < 0.5)
  r <- r[is.finite(r)]
  q_h <- unname(stats::quantile(r, p, type = 7))
  nf <- fit_normal(r)
  z <- stats::qnorm(p)
  tf <- fit_t(r)
  qt_ <- stats::qt(p, tf$nu)
  # E[X | X < q] for a standard t: -dt(q) / p * (nu + q^2) / (nu - 1)
  es_t_std <- -stats::dt(qt_, tf$nu) / p * (tf$nu + qt_^2) / (tf$nu - 1)
  data.frame(
    model = c("historical", "normal", "student_t"),
    VaR = -c(q_h, nf$mu + nf$sigma * z, tf$mu + tf$scale * qt_),
    ES = -c(mean(r[r <= q_h]),
            nf$mu - nf$sigma * stats::dnorm(z) / p,
            tf$mu + tf$scale * es_t_std)
  )
}

#' Kupiec (1995) unconditional coverage test: x exceedances in n days
#' against nominal probability p.
kupiec_test <- function(x, n, p) {
  stopifnot(x >= 0, x <= n, n > 0, p > 0, p < 1)
  xlogy <- function(a, b) ifelse(a == 0, 0, a * log(b))
  phat <- x / n
  ll0 <- xlogy(n - x, 1 - p) + xlogy(x, p)
  ll1 <- xlogy(n - x, 1 - phat) + xlogy(x, phat)
  lr <- -2 * (ll0 - ll1)
  data.frame(exceedances = x, days = n, expected = n * p,
             rate = phat, LR = lr,
             p_value = stats::pchisq(lr, df = 1, lower.tail = FALSE))
}
