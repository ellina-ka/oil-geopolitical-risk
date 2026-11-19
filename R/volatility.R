# Implied vs realised volatility, and mean reversion of implied volatility.

#' Forward-looking realised volatility over the next h trading days,
#' annualised, in the same units as the input returns (percent).
#' rv[t] uses r[t+1], ..., r[t+h]; it is NA where that window is incomplete.
realised_vol_forward <- function(r, h = 21, days_per_year = 252) {
  n <- length(r)
  out <- rep(NA_real_, n)
  r2 <- r^2
  for (t in seq_len(n - h)) {
    w <- r2[(t + 1):(t + h)]
    if (all(is.finite(w))) out[t] <- sqrt(days_per_year / h * sum(w))
  }
  out
}

#' AR(1) on x by OLS in matrix form, beta = (X'X)^{-1} X'y, with
#' Newey-West standard errors. Returns phi, its CI and the implied
#' half-life of shocks, log(0.5) / log(phi), in days.
#' Caveat: OLS phi is biased towards zero in finite samples when phi is
#' near 1 (Kendall 1954), so the half-life is, if anything, understated.
ar1_half_life <- function(x, lag_nw = 10, level = 0.95) {
  x <- x[is.finite(x)]
  y <- x[-1]
  X <- cbind(1, x[-length(x)])
  beta <- solve(crossprod(X), crossprod(X, y))
  fit <- stats::lm(y ~ X[, 2])
  stopifnot(isTRUE(all.equal(unname(stats::coef(fit)), as.vector(beta))))
  V <- sandwich::NeweyWest(fit, lag = lag_nw, prewhite = FALSE)
  phi <- beta[2]
  se <- sqrt(V[2, 2])
  zq <- stats::qnorm(1 - (1 - level) / 2)
  ci <- phi + c(-1, 1) * zq * se
  hl <- function(p) ifelse(p > 0 & p < 1, log(0.5) / log(p), NA_real_)
  data.frame(phi = phi, se_nw = se, phi_lo = ci[1], phi_hi = ci[2],
             half_life = hl(phi), half_life_lo = hl(ci[1]), half_life_hi = hl(ci[2]),
             n = length(y))
}
