# GARCH(1,1) with standardised Student-t innovations, estimated by maximum
# likelihood with optim(). Written from scratch so every assumption is visible:
#   r_t = mu + eps_t,  eps_t = sigma_t * z_t,  z_t ~ t_nu scaled to unit variance
#   sigma2_t = omega + alpha * eps_{t-1}^2 + beta * sigma2_{t-1}
# Constraints: omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1, nu > 2.
# sigma2_1 is set to the sample variance (a common, stated convention).

garch_filter <- function(eps, omega, alpha, beta, sigma2_init = stats::var(eps)) {
  n <- length(eps)
  s2 <- numeric(n)
  s2[1] <- sigma2_init
  for (t in seq_len(n)[-1]) s2[t] <- omega + alpha * eps[t - 1]^2 + beta * s2[t - 1]
  s2 # s2[t] uses information up to t-1 only: an ex-ante forecast
}

# log density of a unit-variance Student-t
std_t_logdens <- function(z, nu) {
  lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * (nu - 2)) -
    (nu + 1) / 2 * log1p(z^2 / (nu - 2))
}

garch_nll <- function(theta, eps, sigma2_init) {
  omega <- theta[1]; alpha <- theta[2]; beta <- theta[3]; nu <- theta[4]
  if (omega <= 0 || alpha < 0 || beta < 0 || alpha + beta >= 1 || nu <= 2) return(1e10)
  s2 <- garch_filter(eps, omega, alpha, beta, sigma2_init)
  -sum(std_t_logdens(eps / sqrt(s2), nu) - 0.5 * log(s2))
}

# unconstrained -> constrained parameters
garch_untransform <- function(p) {
  persistence <- stats::plogis(p[2]); share <- stats::plogis(p[3])
  c(omega = exp(p[1]), alpha = persistence * share,
    beta = persistence * (1 - share), nu = 2 + exp(p[4]))
}

fit_garch_t <- function(r) {
  if (any(!is.finite(r))) stop("fit_garch_t: r must be finite (drop NAs first)")
  if (length(r) < 500) warning("fit_garch_t: fewer than 500 observations; estimates will be noisy")
  mu <- mean(r)
  eps <- r - mu
  s2init <- stats::var(eps)
  p0 <- c(log(0.05 * s2init), stats::qlogis(0.95), stats::qlogis(0.07), log(4))
  obj <- function(p) garch_nll(garch_untransform(p), eps, s2init)
  opt <- stats::optim(p0, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-10))
  if (opt$convergence != 0) warning("fit_garch_t: optim did not converge (code ", opt$convergence, ")")
  theta <- garch_untransform(opt$par)

  # Standard errors from the numerical Hessian on the natural scale.
  se <- tryCatch({
    H <- stats::optimHess(theta, garch_nll, eps = eps, sigma2_init = s2init)
    sqrt(diag(solve(H)))
  }, error = function(e) {
    warning("fit_garch_t: Hessian not invertible, SEs set to NA")
    rep(NA_real_, 4)
  })
  names(se) <- names(theta)
  pers <- unname(theta["alpha"] + theta["beta"])
  list(mu = mu, coef = theta, se = se, loglik = -opt$value,
       sigma2_init = s2init, convergence = opt$convergence,
       sigma = sqrt(garch_filter(eps, theta[["omega"]], theta[["alpha"]],
                                 theta[["beta"]], s2init)),
       persistence = pers,
       half_life_days = log(0.5) / log(pers),
       unconditional_sd = sqrt(theta[["omega"]] / (1 - pers)))
}

#' Ex-ante conditional sd for new data using fixed (already estimated)
#' parameters: used for the out-of-sample VaR backtest.
garch_forecast_sd <- function(r, fit) {
  eps <- r - fit$mu
  sqrt(garch_filter(eps, fit$coef[["omega"]], fit$coef[["alpha"]],
                    fit$coef[["beta"]], fit$sigma2_init))
}

simulate_garch_t <- function(n, omega, alpha, beta, nu, mu = 0, seed = 1) {
  set.seed(seed)
  z <- stats::rt(n, nu) * sqrt((nu - 2) / nu)
  s2 <- numeric(n); eps <- numeric(n)
  s2[1] <- omega / (1 - alpha - beta)
  eps[1] <- sqrt(s2[1]) * z[1]
  for (t in 2:n) {
    s2[t] <- omega + alpha * eps[t - 1]^2 + beta * s2[t - 1]
    eps[t] <- sqrt(s2[t]) * z[t]
  }
  mu + eps
}

# Quantile of the unit-variance t, for VaR.
std_t_quantile <- function(p, nu) stats::qt(p, nu) * sqrt((nu - 2) / nu)
