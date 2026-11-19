# Oil, volatility and Middle East shocks

I work on Gulf political economy, and I wanted to check something I read in market commentary all the time: that oil "prices in" geopolitical risk. This repo looks at two questions:

1. How do Brent returns behave around major Middle East events since 2019?
2. Does oil implied volatility (OVX) overprice or underprice the volatility that follows?

Everything runs from public FRED data in R.

```
Rscript run_all.R             # tables, figures and output/results.md
Rscript run_all.R --refresh   # same, re-downloading the data
Rscript tests/testthat.R      # unit tests
```

## Data

- Brent spot `DCOILBRENTEU` and WTI spot `DCOILWTICO` (EIA via FRED)
- OVX `OVXCLS` (CBOE crude oil volatility index)
- S&P 500 `SP500`, VIX `VIXCLS`, 10y Treasury `DGS10`, broad dollar index `DTWEXBGS`
- `data/events.csv`: 14 events I coded by hand, with a source for each

Raw downloads go into `data/raw/` and are not modified. `MANIFEST.csv` logs when each file was downloaded, its row count and md5 hash. Sample: Jan 2010 to Nov 2026.

## What the code does

**Return distribution.** Moments, a Jarque-Bera test, and normal vs Student-t fits by maximum likelihood. Also 99% VaR and expected shortfall (historical, normal, t).

**GARCH(1,1) with t errors.** I wrote the likelihood myself rather than use rugarch, mainly to understand it. Constraints (alpha + beta < 1, nu > 2) are handled by reparameterising, and standard errors come from the numerical Hessian. One test simulates a GARCH process and checks that the estimator gets the parameters back.

**VaR backtest.** Parameters are estimated on 2010-2019 and then frozen. From 2020 on, each day's VaR uses only information up to the day before. Coverage is checked with the Kupiec test.

**Event study.** Returns are divided by the GARCH conditional volatility before anything else. Otherwise an event that happens in an already volatile week looks bigger than it is. Instead of trusting the N(0,1) reference, I compare the events with 10,000 random sets of ordinary days (placebo test). I also rerun it without events that overlap each other (e.g. 22 and 24 June 2025), using four different windows.

**Implied vs realised.** OVX on day t against WTI realised vol over the next 21 trading days. The event-day gap is compared with random days, and with random days drawn from within a year of each event. That second check is needed because many of the events fall in 2024-25.

**Mean reversion.** An AR(1) on log OVX (OLS, Newey-West SEs), and the half-life it implies.

## Results

Full tables are in `output/results.md`, and figures are in `output/figures/`.

What the data shows directly:
- Brent returns are very fat-tailed. Excess kurtosis is 5.0 excluding Mar-May 2020 and 89 including it. The t distribution (nu about 2.9) fits far better than the normal, which puts 99% expected shortfall at 7.4% against 12.1% in the data.
- Volatility is persistent: alpha + beta = 0.993, a half-life of about 100 trading days.
- Out of sample (2020 to Nov 2026), the GARCH-t 99% VaR was breached 29 times against 17.4 expected (Kupiec p = 0.011). A static normal VaR was breached 93 times. So GARCH-t is much better, but it still under-covers.
- On event days the mean absolute standardised move is 1.75, against a median of 0.72 for random days (placebo p < 0.001). This holds for the 11 non-overlapping events. Over days 0 to 5 the difference is no longer significant (p = 0.09, or 0.16 without overlaps).
- OVX was above subsequent realised volatility on 71% of days, with a median gap of 4.7 points. On event days the median gap was 13.5 points, which is outside the 95% range of both placebo tests.

How I read it:
- The market seems to absorb this kind of news within a day or so. I don't see much drift afterwards.
- The direction depends on what the event means for supply. The biggest falls were on days the market read as de-escalation or a contained response (Oct 2024, and 23-24 June 2025).
- After shocks, options price in more volatility than actually follows. That is consistent with more demand for hedging, or a higher volatility risk premium, after salient events.

What I'm not sure about:
- 14 events is a small sample, and the bootstrap bands are wide.
- I picked the events after the fact. That probably favours events that moved the market. A rule-based list, e.g. from a news-intensity index, would be better.
- Weekend events are mapped to the next trading day. Brent spot is a London price, while OVX and the US series close in New York, so "day 0" is not exactly the same moment for every asset.
- OVX is built on USO options (WTI futures), so it is not a clean measure for Brent.
- EIA Brent spot printed $9.12 on 21 April 2020, far below where Brent futures traded. I report the moments with and without that period.
- 13 May 2019 was also a US-China tariff escalation day, so the cross-asset moves on E01 are not just about the Gulf.
- The full-sample GARCH fit, used to standardise returns, includes the event days themselves.

Things I'd like to add: a rule-based event list, Brent options data instead of OVX, and a Christoffersen independence test for the VaR breaches.

## Files

```
R/           data.R, returns.R, garch.R, events.R, volatility.R
tests/       testthat tests
data/        events.csv, raw/ (FRED downloads + manifest)
output/      results.md, tables/, figures/ (all generated)
run_all.R
```

R 4.2+, with ggplot2, patchwork, sandwich, knitr and testthat.
