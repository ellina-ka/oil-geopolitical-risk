# Reproduces every table and figure from raw FRED data.
# Usage (from the project root):  Rscript run_all.R            # uses cached raw data
#                                  Rscript run_all.R --refresh  # re-downloads
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

refresh <- "--refresh" %in% commandArgs(trailingOnly = TRUE)
START <- as.Date("2010-01-01")
TRAIN_END <- as.Date("2019-12-31") # VaR model estimated up to here, tested after
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
save_tab <- function(x, name) utils::write.csv(x, file.path("output/tables", paste0(name, ".csv")), row.names = FALSE)

# ---- 1. Data --------------------------------------------------------------
ids <- c(brent = "DCOILBRENTEU", wti = "DCOILWTICO", ovx = "OVXCLS",
         vix = "VIXCLS", spx = "SP500", ust10 = "DGS10", usd = "DTWEXBGS")
raw <- lapply(ids, fetch_fred, refresh = refresh)
events <- read_events()

brent <- subset(raw$brent, date >= START)
brent$r <- log_returns(brent$value)
brent <- brent[is.finite(brent$r), ]

# ---- 2. Return distribution ----------------------------------------------
desc <- describe_returns(brent$r)
nf <- fit_normal(brent$r); tf <- fit_t(brent$r)
dist_fit <- data.frame(model = c("normal", "student_t"),
                       loglik = c(nf$loglik, tf$loglik), aic = c(nf$aic, tf$aic),
                       nu = c(NA, tf$nu))
tails <- var_es(brent$r, p = 0.01)
covid <- brent$date >= as.Date("2020-03-01") & brent$date <= as.Date("2020-05-31")
desc <- rbind(cbind(sample = "full", desc),
              cbind(sample = "excl. Mar-May 2020", describe_returns(brent$r[!covid])))
save_tab(desc, "01_brent_moments"); save_tab(dist_fit, "02_distribution_fit"); save_tab(tails, "03_var_es_99")

# ---- 3. GARCH(1,1)-t and out-of-sample VaR backtest -----------------------
g <- fit_garch_t(brent$r)
brent$sigma <- g$sigma
brent$z <- (brent$r - g$mu) / brent$sigma
garch_tab <- data.frame(parameter = names(g$coef), estimate = unname(g$coef), se = unname(g$se))
garch_tab <- rbind(garch_tab, data.frame(
  parameter = c("persistence", "half_life_days", "unconditional_annual_vol_pct"),
  estimate = c(g$persistence, g$half_life_days, g$unconditional_sd * sqrt(252)), se = NA))
save_tab(garch_tab, "04_garch_t_full_sample")

train <- brent$date <= TRAIN_END
g_train <- fit_garch_t(brent$r[train])
sd_ex_ante <- garch_forecast_sd(brent$r, g_train) # parameters fixed at 2019 values
test <- !train
p <- 0.01
var_garch <- g_train$mu + sd_ex_ante * std_t_quantile(p, g_train$coef[["nu"]])
nf_train <- fit_normal(brent$r[train])
var_static <- nf_train$mu + nf_train$sigma * stats::qnorm(p)
backtest <- rbind(
  cbind(model = "GARCH(1,1)-t, 2010-19 parameters", kupiec_test(sum(brent$r[test] < var_garch[test]), sum(test), p)),
  cbind(model = "Static normal, 2010-19 parameters", kupiec_test(sum(brent$r[test] < var_static), sum(test), p)))
save_tab(backtest, "05_var_backtest_2020_onwards")

# ---- 4. Event study -------------------------------------------------------
windows <- list(d0 = c(0, 0), d0_1 = c(0, 1), d0_5 = c(0, 5), dm1_1 = c(-1, 1))
ev_tab <- event_table(brent, events, windows)
save_tab(ev_tab, "06_event_returns")

placebo <- do.call(rbind, lapply(windows, function(w) placebo_test(brent, events, w[1], w[2])))
indep <- events[!ev_tab$overlaps_previous, ]
placebo_indep <- do.call(rbind, lapply(windows, function(w) placebo_test(brent, indep, w[1], w[2], exclude_dates = events$event_date)))
placebo$sample <- "all events"; placebo_indep$sample <- "non-overlapping events only"
WAR_START <- as.Date("2026-02-28")
by_period <- lapply(list("2019-2025 events" = events$event_date < WAR_START,
                         "2026 war events" = events$event_date >= WAR_START), function(keep) {
  do.call(rbind, lapply(windows, function(w) placebo_test(brent, events[keep, ], w[1], w[2], exclude_dates = events$event_date)))
})
for (nm in names(by_period)) by_period[[nm]]$sample <- nm
war_within <- do.call(rbind, lapply(windows, function(w)
  placebo_test(brent, events[events$event_date >= WAR_START, ], w[1], w[2],
               exclude_dates = events$event_date, pool_from = WAR_START)))
war_within$sample <- "2026 war events vs other war days"
war_within5 <- do.call(rbind, lapply(windows, function(w)
  placebo_test(brent, events[events$event_date >= WAR_START, ], w[1], w[2], exclusion = 3,
               exclude_dates = events$event_date, pool_from = WAR_START)))
war_within5$sample <- "2026 war events vs other war days (+/-3 day exclusion)"
placebo_all <- rbind(placebo, placebo_indep, by_period[[1]], by_period[[2]], war_within, war_within5)

# Before vs during the 2026 war: levels of risk, not event effects.
regime <- function(keep_b, label) data.frame(
  period = label,
  days = sum(keep_b),
  brent_ann_vol_pct = sd(brent$r[keep_b]) * sqrt(252),
  mean_abs_daily_move_pct = mean(abs(brent$r[keep_b])),
  share_days_move_over_5pct = mean(abs(brent$r[keep_b]) > 5),
  mean_garch_vol_pct = mean(brent$sigma[keep_b]) * sqrt(252))
regime_tab <- rbind(
  regime(brent$date >= as.Date("2019-01-01") & brent$date < WAR_START, "Jan 2019 to Feb 2026"),
  regime(brent$date >= WAR_START, "Since 28 Feb 2026 (war)"))
save_tab(regime_tab, "07b_war_regime")
save_tab(placebo_all, "07_placebo_tests")

# ---- 5. Cross-asset reaction on day 0 -------------------------------------
lvl <- function(id) raw[[id]][, c("date", "value")]
chg <- function(d, type) {
  d$chg <- if (type == "log") log_returns(d$value) else c(NA, diff(d$value))
  d[, c("date", "chg")]
}
cross <- Reduce(function(a, b) merge(a, b, by = "date", all.x = TRUE), list(
  brent[, c("date", "r")],
  setNames(chg(lvl("spx"), "log"), c("date", "spx_pct")),
  setNames(chg(lvl("usd"), "log"), c("date", "usd_pct")),
  setNames(chg(lvl("vix"), "diff"), c("date", "vix_pts")),
  setNames(transform(chg(lvl("ust10"), "diff"), chg = 100 * chg), c("date", "ust10_bp")),
  setNames(chg(lvl("ovx"), "diff"), c("date", "ovx_pts"))))
names(cross)[2] <- "brent_pct"
cross_ev <- merge(ev_tab[, c("event_id", "day0", "direction")], cross,
                  by.x = "day0", by.y = "date", all.x = TRUE)
cross_ev <- cross_ev[order(cross_ev$day0), ]
save_tab(cross_ev, "08_cross_asset_day0")

# ---- 6. Implied vs realised volatility (WTI, because OVX is on WTI-linked options)
wti <- subset(raw$wti, date >= START)
wti$r <- log_returns(wti$value, on_nonpositive = "na")
vol <- merge(wti[, c("date", "r")], setNames(lvl("ovx"), c("date", "ovx")), by = "date")
vol$rv21 <- realised_vol_forward(vol$r, h = 21)
vol$gap <- vol$ovx - vol$rv21 # implied minus subsequently realised, vol points
vrp <- data.frame(
  n_days = sum(is.finite(vol$gap)),
  mean_ovx = mean(vol$ovx[is.finite(vol$gap)]),
  mean_rv21 = mean(vol$rv21, na.rm = TRUE),
  mean_gap = mean(vol$gap, na.rm = TRUE),
  median_gap = stats::median(vol$gap, na.rm = TRUE),
  share_implied_above_realised = mean(vol$gap > 0, na.rm = TRUE))
save_tab(vrp, "09_implied_vs_realised")

i0_vol <- map_day0(events$event_date, vol$date)
ev_vol <- data.frame(event_id = events$event_id, day0 = vol$date[i0_vol],
                     ovx_before = vol$ovx[i0_vol - 1], ovx_day0 = vol$ovx[i0_vol],
                     ovx_day5 = vol$ovx[i0_vol + 5], ovx_day20 = vol$ovx[i0_vol + 20],
                     rv21_after_day0 = vol$rv21[i0_vol])
ev_vol$gap_day0 <- ev_vol$ovx_day0 - ev_vol$rv21_after_day0
save_tab(ev_vol, "10_ovx_around_events")

# Placebo: is the implied-minus-realised gap on event days unusual?
# Two comparisons per period: any non-event day, and a day matched to each
# event within +/-250 trading days (so events clustering in high- or
# low-premium years do not drive the result). Two-sided p-values, because
# the sign of the gap differs before and during the 2026 war.
ok_gap <- which(is.finite(vol$gap))
near <- unique(unlist(lapply(i0_vol[!is.na(i0_vol)], function(i) (i - 10):(i + 10))))
pool_gap <- setdiff(ok_gap, near)
gap_placebo <- function(sel, label, n_draws = 10000, seed = 11) {
  keep <- sel & is.finite(ev_vol$gap_day0)
  obs <- stats::median(ev_vol$gap_day0[keep])
  set.seed(seed)
  any_day <- replicate(n_draws, stats::median(vol$gap[sample(pool_gap, sum(keep))]))
  local_pool <- lapply(i0_vol[keep], function(i) intersect(pool_gap, (i - 250):(i + 250)))
  matched <- replicate(n_draws, stats::median(vapply(local_pool, function(p)
    vol$gap[p[sample.int(length(p), 1)]], numeric(1))))
  two_sided <- function(d) min(1, 2 * min((1 + sum(d >= obs)) / (n_draws + 1),
                                          (1 + sum(d <= obs)) / (n_draws + 1)))
  data.frame(events = label, placebo = c("any non-event day", "matched within 250 days"),
             n_events = sum(keep), observed_median_gap = obs,
             placebo_median = c(stats::median(any_day), stats::median(matched)),
             placebo_p05 = c(stats::quantile(any_day, 0.05), stats::quantile(matched, 0.05)),
             placebo_p95 = c(stats::quantile(any_day, 0.95), stats::quantile(matched, 0.95)),
             p_value_two_sided = c(two_sided(any_day), two_sided(matched)))
}
gap_test <- rbind(gap_placebo(rep(TRUE, nrow(ev_vol)), "all"),
                  gap_placebo(events$event_date < WAR_START, "2019-2025"),
                  gap_placebo(events$event_date >= WAR_START, "2026 war"))
rownames(gap_test) <- NULL
save_tab(gap_test, "10b_gap_placebo")

ovx_ar <- ar1_half_life(log(raw$ovx$value[raw$ovx$date >= START]))
save_tab(ovx_ar, "11_ovx_ar1_half_life")

# ---- 7. Figures -----------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
ev_day0 <- ev_tab$day0

zz <- (brent$r - mean(brent$r)) / sd(brent$r)
qq <- data.frame(theoretical = stats::qnorm(stats::ppoints(length(zz))), sample = sort(zz))
f1 <- ggplot(qq, aes(theoretical, sample)) + geom_abline(colour = "grey60") +
  geom_point(size = 0.6, alpha = 0.5) +
  labs(title = "Brent daily returns vs normal", x = "Normal quantile", y = "Standardised return")
f1b <- ggplot(brent, aes(r)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, fill = "grey75") +
  stat_function(fun = function(x) stats::dnorm(x, nf$mu, nf$sigma), aes(colour = "Normal"), n = 1000) +
  stat_function(fun = function(x) stats::dt((x - tf$mu) / tf$scale, tf$nu) / tf$scale, aes(colour = "Student-t"), n = 1000) +
  coord_cartesian(xlim = c(-8, 8)) + scale_colour_manual(values = c("#d95f02", "#1b9e77")) +
  labs(title = sprintf("Fitted densities (t: nu = %.1f)", tf$nu), x = "Daily log return (%)", y = NULL, colour = NULL)
ggsave("output/figures/01_return_distribution.png", f1 + f1b, width = 10, height = 4, dpi = 150)

f2 <- ggplot(brent, aes(date, sigma * sqrt(252))) + geom_line(linewidth = 0.3) +
  geom_vline(xintercept = ev_day0, colour = "#d95f02", alpha = 0.5, linetype = 2) +
  scale_y_log10() +
  labs(title = "Brent GARCH(1,1)-t conditional volatility (annualised, log scale)",
       subtitle = "Dashed lines: geopolitical events in data/events.csv", x = NULL, y = "%")
ggsave("output/figures/02_garch_volatility.png", f2, width = 10, height = 4, dpi = 150)

i0_b <- map_day0(events$event_date, brent$date)
path_z <- event_path(cumsum(ifelse(is.finite(brent$z), brent$z, 0)), i0_b, -5, 10)
path_z$mean <- path_z$mean - path_z$mean[path_z$k == -1]
path_z$lo <- path_z$lo - path_z$lo[path_z$k == -1]; path_z$hi <- path_z$hi - path_z$hi[path_z$k == -1]
abs_z <- event_path(abs(brent$z), i0_b, -5, 10)
f3 <- ggplot(abs_z, aes(k, mean)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_line() + geom_point() + geom_hline(yintercept = sqrt(2 / pi) * 1, linetype = 3) +
  labs(title = "|Standardised abnormal return| around events",
       subtitle = "Band: 90% bootstrap over events.\nDotted: expected |z| under N(0,1) noise.",
       x = "Trading days from event", y = "|z|")
ovx_path <- event_path(log(vol$ovx), i0_vol, -5, 30)
base <- ovx_path$mean[ovx_path$k == -1]
f4 <- ggplot(ovx_path, aes(k, 100 * (mean - base))) +
  geom_ribbon(aes(ymin = 100 * (lo - base), ymax = 100 * (hi - base)), fill = "grey85") +
  geom_line() + geom_hline(yintercept = 0, linetype = 3) +
  labs(title = "OVX around events (% vs day -1)", subtitle = sprintf("Band: 90%% bootstrap over events.\nUnconditional AR(1) half-life of log OVX: %.0f trading days", ovx_ar$half_life),
       x = "Trading days from event", y = "% change")
ggsave("output/figures/03_event_paths.png", f3 + f4, width = 10, height = 4, dpi = 150)

f5 <- ggplot(subset(vol, is.finite(rv21)), aes(date)) +
  geom_line(aes(y = ovx, colour = "OVX (implied, 30d)"), linewidth = 0.3) +
  geom_line(aes(y = rv21, colour = "WTI realised, next 21 days"), linewidth = 0.3) +
  scale_y_log10() + scale_colour_manual(values = c("#1b9e77", "#7570b3")) +
  labs(title = "Implied vs subsequently realised oil volatility", x = NULL, y = "% annualised (log scale)", colour = NULL) +
  theme(legend.position = "bottom")
ggsave("output/figures/04_implied_vs_realised.png", f5, width = 10, height = 4, dpi = 150)

# ---- 8. Results file ------------------------------------------------------
kab <- function(x, digits = 3) paste(knitr::kable(x, digits = digits, row.names = FALSE), collapse = "\n")
md <- c("# Results (generated by run_all.R)", "",
        sprintf("Brent sample: %s to %s (%d daily returns). Generated %s.",
                min(brent$date), max(brent$date), nrow(brent), Sys.Date()), "",
        "## Return distribution", kab(desc), "", kab(dist_fit), "",
        "### 1-day 99% VaR and ES (loss, %)", kab(tails), "",
        "## GARCH(1,1)-t (full sample)", kab(garch_tab, 4), "",
        "## Out-of-sample 99% VaR backtest (2020 onwards; parameters frozen at 2019)", kab(backtest), "",
        "## Event returns (% log return; SCAR = standardised cumulative abnormal return)", kab(ev_tab, 2), "",
        "## Placebo tests (mean |SCAR| vs random non-event days)", kab(placebo_all), "",
        "## Before vs during the 2026 Iran war", kab(regime_tab), "",
        "## Cross-asset moves on day 0", kab(cross_ev, 2), "",
        "## Implied vs realised (WTI / OVX)", kab(vrp, 2), "", kab(ev_vol, 1), "", "Placebo for the day-0 gap:", kab(gap_test, 3), "",
        "## Mean reversion of log OVX (AR(1), Newey-West SE)", kab(ovx_ar, 4), "",
        "## Session", "```", utils::capture.output(print(utils::sessionInfo())), "```")
writeLines(md, "output/results.md")
message("Done: output/results.md, output/tables/, output/figures/")
