# Event study on GARCH-standardised returns.
#
# Abnormal return z_t = (r_t - mu) / sigma_t, with sigma_t the ex-ante GARCH
# conditional sd. Scaling by sigma_t (not a constant sd) stops events that
# land in already-volatile periods from looking more significant than they are.
# Window statistic: SCAR[a,b] = sum(z_t, t = a..b) / sqrt(b - a + 1), which is
# approximately N(0,1) under "no event effect" if the GARCH model is adequate.
# Because that "if" is a strong assumption, inference also uses a placebo
# distribution built from non-event days.

#' Index of the first trading date on or after each event date.
#' Weekend and holiday events therefore map to the next trading day.
map_day0 <- function(event_dates, trading_dates) {
  stopifnot(!is.unsorted(trading_dates))
  idx <- findInterval(as.numeric(event_dates) - 1, as.numeric(trading_dates)) + 1
  idx[idx > length(trading_dates)] <- NA
  idx
}

window_sum <- function(x, i0, a, b) {
  vapply(i0, function(i) {
    lo <- i + a; hi <- i + b
    if (is.na(i) || lo < 1 || hi > length(x)) return(NA_real_)
    sum(x[lo:hi])
  }, numeric(1))
}

#' Flag events whose day 0 falls within `gap` trading days of an earlier
#' event: their windows overlap, so they are not independent observations.
flag_overlap <- function(i0, gap = 5) {
  o <- order(i0)
  prev <- c(-Inf, i0[o][-length(o)])
  flag <- logical(length(i0))
  flag[o] <- (i0[o] - prev) <= gap
  flag
}

event_table <- function(df, events, windows) {
  # df: date, r (pct log return), z (standardised abnormal return)
  i0 <- map_day0(events$event_date, df$date)
  out <- data.frame(event_id = events$event_id,
                    event_date = events$event_date,
                    day0 = df$date[i0],
                    direction = events$direction,
                    overlaps_previous = flag_overlap(i0))
  for (w in names(windows)) {
    a <- windows[[w]][1]; b <- windows[[w]][2]
    out[[paste0("ret_", w)]] <- window_sum(df$r, i0, a, b)
    out[[paste0("scar_", w)]] <- window_sum(df$z, i0, a, b) / sqrt(b - a + 1)
  }
  out
}

#' Placebo test: is the mean |SCAR| across events larger than for the same
#' number of randomly drawn non-event days? Returns the share of placebo
#' draws at least as large as the observed statistic (one-sided p-value).
placebo_test <- function(df, events, a, b, n_draws = 10000, exclusion = 10, seed = 42,
                         exclude_dates = events$event_date, pool_from = NULL) {
  i0 <- map_day0(events$event_date, df$date)
  L <- b - a + 1
  scar_all <- window_sum(df$z, seq_len(nrow(df)), a, b) / sqrt(L)
  # Events too recent for the full window (e.g. [0,5] at the end of the
  # sample) are dropped from this window only, so the placebo uses the same N.
  i0 <- i0[!is.na(i0)]
  i0 <- i0[is.finite(scar_all[i0])]
  # Placebo days stay away from every listed event, not just the ones tested.
  ix <- map_day0(exclude_dates, df$date)
  ix <- ix[!is.na(ix)]
  near_event <- unique(unlist(lapply(ix, function(i) (i - exclusion):(i + exclusion))))
  pool <- setdiff(which(is.finite(scar_all)), near_event)
  # Optionally draw placebo days only from a sub-period (e.g. inside the war),
  # so events are compared with ordinary days of the same regime.
  if (!is.null(pool_from)) pool <- pool[df$date[pool] >= pool_from]
  if (length(pool) < 5 * length(i0)) warning("placebo pool is small: ", length(pool), " days")
  observed <- mean(abs(scar_all[i0]), na.rm = TRUE)
  set.seed(seed)
  draws <- replicate(n_draws, mean(abs(scar_all[sample(pool, length(i0))])))
  data.frame(window = sprintf("[%d,%d]", a, b), n_events = length(i0),
             observed_mean_abs_scar = observed,
             placebo_median = stats::median(draws),
             placebo_p95 = unname(stats::quantile(draws, 0.95)),
             p_value = (1 + sum(draws >= observed)) / (n_draws + 1))
}

#' Mean event-time path of a series with a percentile bootstrap band
#' (resampling events). With ~14 events the band is wide by construction.
event_path <- function(x, i0, from, to, n_boot = 5000, level = 0.90, seed = 7) {
  i0 <- i0[!is.na(i0) & i0 + from >= 1 & i0 + to <= length(x)]
  k <- from:to
  M <- t(vapply(i0, function(i) x[i + k], numeric(length(k))))
  set.seed(seed)
  boot <- replicate(n_boot, colMeans(M[sample(nrow(M), replace = TRUE), , drop = FALSE]))
  alpha <- (1 - level) / 2
  data.frame(k = k, mean = colMeans(M),
             lo = apply(boot, 1, stats::quantile, alpha),
             hi = apply(boot, 1, stats::quantile, 1 - alpha),
             n_events = nrow(M))
}
