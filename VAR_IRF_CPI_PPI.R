# CBF CPI/PPI scenarios
# Created by Robert Bishop
# Very decent amount of use of Codex

settings = list(
  holdings_path =  "C:/Users/Roberto/Downloads/RJFPortfolio 9.26.2026 (1).csv", # Whatever path on ur computer for the RJ CSV holdings template just link it here
  holdings_as_of = "2026-09-24", # Ex. "2026-09-24"
  start_date = "2010-01-01",
  estimation_end_month = "2025-09",
  target_month = "2026-10",
  cpi_mom_pct = 0.246939616833026,                 # Ex. assuntion
  ppi_mom_pct = 0.267626761776557,
  seed = 20261006L,
  irf_runs = 1000L,
  backtest_months = 24L,
  output_dir = "outputs/var_irf_test"
)

forecast_reference = data.frame(
  month = "2026-10",
  CPI = 0.246939616833026,
  PPI = 0.267626761776557
)

month_date = function(x) as.Date(paste0(x, "-01"))

to_log_pct = function(x) {
  if (any(!is.finite(x)) || any(x <= -100))
    stop("Changes must be finite and greater than -100%.")
  100 * log1p(x / 100)
}

to_simple_pct = function(x) 100 * expm1(x / 100)

check_months = function(dates) {
  if (!length(dates) || anyNA(dates) || anyDuplicated(dates) ||
      !identical(as.character(dates),
                 as.character(seq(min(dates), max(dates), by = "month"))))
    stop("History has duplicate, missing, or unordered calendar months.")
}

monthly_changes = function(dates, levels, name) {
  dates = as.Date(format(as.Date(dates), "%Y-%m-01"))
  check_months(dates)
  invalid = !is.finite(levels) | levels <= 0
  if (any(invalid))
    stop(paste(name, "invalid prices:", paste(format(dates[invalid], "%Y-%m"), collapse = ", ")))
  result = data.frame(month = dates[-1L], value = 100 * diff(log(levels)))
  names(result)[2L] = name
  result
}

read_holdings = function(path = NULL, as_of = NULL) {
  if (is.null(path)) path = file.choose()
  if (is.null(as_of)) as_of = readline("Holdings valuation date (YYYY-MM-DD): ")
  as_of = as.character(as.Date(as_of))
  if (length(as_of) != 1L || is.na(as_of)) stop("Enter the holdings valuation date.")
  encoding = if (all(validUTF8(readLines(path, warn = FALSE)))) "UTF-8" else "latin1"
  csv = read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, encoding = encoding)
  names(csv) = sub("^\ufeff", "", names(csv), useBytes = TRUE)
  columns = c("SYMBOL/CUSIP", "Current Value", "Product Type")
  if (!all(columns %in% names(csv)) || !nrow(csv)) stop("Use a populated Raymond James CSV.")
  types = c("Stock" = "stock", "Funds" = "fund", "Cash & Cash Alternatives" = "cash")
  p = data.frame(
    symbol = toupper(trimws(csv[["SYMBOL/CUSIP"]])),
    assetType = unname(types[trimws(csv[["Product Type"]])]),
    value_usd = as.numeric(gsub("[,$[:space:]]", "", csv[["Current Value"]]))
  )
  if (anyNA(p$assetType) || any(!is.finite(p$value_usd)) ||
      any(p$value_usd < 0) || sum(p$value_usd) <= 0) stop("Check product types and current values.")
  securities = p$assetType != "cash"
  if (anyNA(p$symbol[securities]) || any(!nzchar(p$symbol[securities])) ||
      any(grepl("^[A-Z0-9]{9}$", p$symbol[securities]))) stop("Use tickers rather than CUSIPs.")
  if (anyDuplicated(p$symbol[securities])) stop("Consolidate duplicate tickers before importing.")
  list(as_of = as_of, positions = p, path = normalizePath(path, winslash = "/"))
}

fetch_macro = function(start_date, end_month = NULL) {
  fetch = function(symbol, name) {
    x = quantmod::getSymbols(symbol, src = "FRED", auto.assign = FALSE)
    first = seq(as.Date(format(as.Date(start_date), "%Y-%m-01")),
                length.out = 2L, by = "-1 month")[2L]
    last = if (is.null(end_month))
      as.Date(format(as.Date(Sys.time(), tz = "America/New_York"), "%Y-%m-01")) - 1L else
        seq(month_date(end_month), length.out = 2L, by = "month")[2L] - 1L
    x = x[paste0(first, "/", last)]
    monthly_changes(zoo::index(x), as.numeric(x), name)
  }
  cpi = fetch("CPIAUCSL", "CPI")
  ppi = fetch("PPIFIS", "PPI")
  result = merge(ppi, cpi, by = "month", sort = TRUE)
  result = result[result$month >= as.Date(start_date) &
                    result$month < as.Date(format(as.Date(Sys.time(), tz = "America/New_York"), "%Y-%m-01")), ]
  result
}

fetch_returns = function(symbol, start_date, cutoff) {
  first = seq(as.Date(format(as.Date(start_date), "%Y-%m-01")),
              length.out = 2L, by = "-1 month")[2L]
  last = seq(cutoff, length.out = 2L, by = "month")[2L]
  x = quantmod::getSymbols(symbol, src = "yahoo", from = first,
                          to = last, auto.assign = FALSE)
  adjusted = quantmod::Ad(x)
  dates = as.Date(zoo::index(adjusted))
  months = as.Date(format(dates, "%Y-%m-01"))
  ends = !duplicated(months, fromLast = TRUE)
  first_complete = first
  if (min(dates) > first + 7) {
    first_complete = seq(min(months), length.out = 2L, by = "month")[2L]
  }
  keep = ends & months >= first_complete & months <= cutoff
  if (max(dates) < last - 7L) stop("Yahoo's last month is incomplete at the macro cutoff.")
  monthly_changes(months[keep], as.numeric(adjusted)[keep], "Return")
}

aligned_history = function(macro, returns) {
  result = merge(macro, returns, by = "month", sort = TRUE)
  check_months(result$month)
  if (max(result$month) != max(macro$month))
    stop("Yahoo history does not reach the macro cutoff.")
  result
}

fit_var = function(y) {
  y = as.matrix(y)
  n = nrow(y)
  if (n < 60L) stop(paste("Insufficient history:", n, "months; need 60."))
  eligible = (1L:6L)[n - (1L:6L) >= 10 * (3 * (1L:6L) + 1)]
  selected = vars::VARselect(y, lag.max = max(eligible), type = "const")
  fits = lapply(eligible, function(p) vars::VAR(y, p = p, type = "const"))
  checks = vapply(fits, function(fit) {
    tryCatch({ model_checks(fit); "passed" }, error = function(e) conditionMessage(e))
  }, character(1))
  bic = as.numeric(selected$criteria["SC(n)", ])
  passing = which(checks == "passed")
  if (!length(passing)) stop(paste("No eligible lag passes model checks:", checks[which.min(bic)]))
  fit = fits[[passing[which.min(bic[passing])]]]
  fit$lag_selection = data.frame(lag = eligible, BIC = bic, check = checks)
  fit$call$p = fit$p
  fit
}

residual_covariance = function(model) {
  e = as.matrix(residuals(model))
  df = nrow(e) - (model$K * model$p + 1L)
  sigma = crossprod(e) / df
  if (any(!is.finite(sigma)) || rcond(sigma) < 1e-10)
    stop("Residual covariance is singular or poorly conditioned.")
  chol(sigma)
  sigma
}

model_checks = function(model) {
  if (any(vars::roots(model) >= 1)) stop("VAR is unstable.")
  residual_covariance(model)
  serial = vars::serial.test(model, lags.pt = 12, type = "PT.adjusted")
  pv = as.numeric(serial$serial$p.value)
  if (pv < 0.05)
    stop(sprintf("Residual autocorrelation: adjusted Portmanteau p=%.4f.", pv))
  pv
}

stationarity_notes = function(y) {
  notes = character()
  for (name in colnames(y)) {
    a = suppressWarnings(tseries::adf.test(y[, name])$p.value)
    k = suppressWarnings(tseries::kpss.test(y[, name], null = "Level")$p.value)
    if (a >= 0.05 || k < 0.05)
      notes = c(notes, sprintf("%s stationarity review (ADF %.3f, KPSS %.3f)", name, a, k))
  }
  paste(notes, collapse = "; ")
}

forecast_distribution = function(model, horizon) {
  k = model$K
  sigma = residual_covariance(model)
  phi = vars::Phi(model, nstep = max(1L, horizon - 1L))
  predictions = predict(model, n.ahead = horizon)$fcst
  means = sapply(colnames(model$y), function(name) predictions[[name]][, "fcst"])
  propagation = matrix(0, k * horizon, k * horizon)
  for (h in seq_len(horizon)) {
    rows = ((h - 1L) * k + 1L):(h * k)
    for (j in seq_len(h)) {
      cols = ((j - 1L) * k + 1L):(j * k)
      propagation[rows, cols] = phi[, , h - j + 1L]
    }
  }
  omega = propagation %*% kronecker(diag(horizon), sigma) %*% t(propagation)
  list(mean = as.vector(t(means)), covariance = omega)
}

condition_forecast = function(distribution, indices, values) {
  mu = distribution$mean
  omega = distribution$covariance
  selected = omega[indices, indices, drop = FALSE]
  gain = t(solve(selected, omega[indices, , drop = FALSE]))
  conditional_mean = as.vector(mu + gain %*% (values - mu[indices]))
  covariance = omega - gain %*% omega[indices, , drop = FALSE]
  covariance = (covariance + t(covariance)) / 2
  list(mean = conditional_mean, covariance = covariance)
}

arithmetic_forecast = function(mean, variance) {
  if (!is.finite(mean) || !is.finite(variance) || variance < -1e-7)
    stop("Invalid conditional mean or variance.")
  variance = max(variance, 0)
  result = c(mean = 100 * expm1(mean / 100 + variance / 20000),
    lower = to_simple_pct(mean + qnorm(0.025) * sqrt(variance)),
    upper = to_simple_pct(mean + qnorm(0.975) * sqrt(variance)))
  if (any(!is.finite(result))) stop("Scenario returns exceed the numerical range.")
  result
}

scenario_forecast = function(model, horizon, cpi, ppi) {
  distribution = forecast_distribution(model, horizon)
  indices = (horizon - 1L) * 3L + c(1L, 2L)
  baseline = condition_forecast(distribution, indices, distribution$mean[indices])
  scenario = condition_forecast(distribution, indices, to_log_pct(c(ppi, cpi)))
  r = horizon * 3L
  base_return = arithmetic_forecast(baseline$mean[r], baseline$covariance[r, r])
  scenario_return = arithmetic_forecast(scenario$mean[r], scenario$covariance[r, r])
  c(baseline_ppi_pct = to_simple_pct(distribution$mean[indices[1L]]),
    baseline_cpi_pct = to_simple_pct(distribution$mean[indices[2L]]),
    baseline_return_pct = unname(base_return["mean"]),
    scenario_return_pct = unname(scenario_return["mean"]),
    difference_pp = unname(scenario_return["mean"] - base_return["mean"]),
    lower_95_pct = unname(scenario_return["lower"]),
    upper_95_pct = unname(scenario_return["upper"]))
}

backtest = function(y, months = 24L) {
  n = nrow(y)
  first = max(61L, n - months + 1L)
  rows = list()
  if (first > n) return(data.frame())
  for (j in first:n) {
    training = y[seq_len(j - 1L), , drop = FALSE]
    result = tryCatch({
      fit = fit_var(training)
      model_checks(fit)
      distribution = forecast_distribution(fit, 1L)
      forecast = arithmetic_forecast(distribution$mean[3L], distribution$covariance[3L, 3L])["mean"]
      conditional = scenario_forecast(fit, 1L, to_simple_pct(y[j, "CPI"]),
                                     to_simple_pct(y[j, "PPI"]))["scenario_return_pct"]
      data.frame(index = j, actual_pct = to_simple_pct(y[j, "Return"]),
                 var_pct = unname(forecast),
                 historical_mean_pct = mean(to_simple_pct(training[, "Return"])),
                 ex_post_conditional_pct = unname(conditional), status = "modeled")
    }, error = function(e) data.frame(index = j, actual_pct = to_simple_pct(y[j, "Return"]),
                                      var_pct = NA_real_, historical_mean_pct = NA_real_,
                                      ex_post_conditional_pct = NA_real_, status = conditionMessage(e)))
    rows[[length(rows) + 1L]] = result
  }
  do.call(rbind, rows)
}

holding_irfs = function(model, runs, seed) {
  make = function(fit) {
    fit$call$p = fit$p
    vars::irf(fit, impulse = c("PPI", "CPI"), response = "Return",
                                n.ahead = 12, ortho = TRUE, cumulative = FALSE,
                                boot = TRUE, ci = 0.95, runs = runs, seed = seed)
  }
  reverse = vars::VAR(model$y[, c("CPI", "PPI", "Return")], p = model$p, type = "const")
  list(PPI_first = make(model), CPI_first = make(reverse))
}

plot_holding_irfs = function(irfs, ticker, cutoff = NULL) {
  old = par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  on.exit(par(old))
  for (ordering in names(irfs)) {
    x = irfs[[ordering]]
    for (impulse in c("PPI", "CPI")) {
      center = as.numeric(x$irf[[impulse]])
      lower = as.numeric(x$Lower[[impulse]])
      upper = as.numeric(x$Upper[[impulse]])
      h = 0:12
      plot(h, center, type = "n", ylim = range(c(lower, upper, 0)),
           xlab = "Months after innovation", ylab = "Adjusted-price log return (pp)",
           main = paste(impulse, "innovation |", gsub("_", " ", ordering)), cex.main = 0.95)
      polygon(c(h, rev(h)), c(lower, rev(upper)), col = "grey90", border = NA)
      abline(h = 0, lty = 3)
      lines(h, center, col = "steelblue4", lwd = 2)
    }
  }
  label = paste(ticker, "| One-standard-deviation innovations | 95% pointwise bootstrap bands")
  if (!is.null(cutoff)) label = paste(label, "| Cutoff", format(as.Date(cutoff), "%Y-%m"))
  mtext(label, outer = TRUE, line = 1, cex = 0.9)
}

plot_scenario_impacts = function(results, coverage, config, cutoff, as_of, kind = "dollar") {
  d = results[results$status == "modeled", ]
  d = d[order(d$dollar_difference), ]
  old = par(mar = c(5, 6, 5, 2), oma = c(4, 0, 0, 0), las = 1)
  on.exit(par(old))
  heading = if (kind == "dollar") "Modeled dollar impact by holding" else "Baseline and scenario returns"
  if (!nrow(d)) {
    plot.new()
    text(0.5, 0.5, "No holdings passed the model checks.")
  } else if (kind == "dollar") {
    span = max(diff(range(c(0, d$dollar_difference))), 1)
    limits = range(c(0, d$dollar_difference)) + c(-1, 1) * span * 0.25
    colors = ifelse(d$dollar_difference >= 0, "#28658C", "#C17838")
    y = barplot(d$dollar_difference, names.arg = d$ticker, horiz = TRUE,
                col = colors, border = NA, xlim = limits,
                xlab = "Scenario minus baseline (USD)", cex.names = 0.9)
    abline(v = 0, col = "grey35")
    labels = paste0(ifelse(d$dollar_difference < 0, "-$", "+$"),
                    formatC(abs(d$dollar_difference), format = "f", digits = 2))
    text(d$dollar_difference, y, labels,
         pos = ifelse(d$dollar_difference < 0, 2, 4), cex = 0.85)
  } else {
    y = seq_len(nrow(d))
    limits = range(c(0, d$baseline_return_pct, d$lower_95_pct, d$upper_95_pct))
    plot(d$scenario_return_pct, y, type = "n", xlim = limits,
         ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", ylab = "",
         xlab = "Target-month adjusted-price return (%)")
    axis(2, at = y, labels = d$ticker, tick = FALSE, cex.axis = 0.9)
    abline(v = 0, col = "grey70", lty = 3)
    segments(d$lower_95_pct, y, d$upper_95_pct, y, col = "#8BA8BA", lwd = 2)
    points(d$baseline_return_pct, y, pch = 1, col = "grey20", cex = 1.1)
    points(d$scenario_return_pct, y, pch = 16, col = "#28658C", cex = 0.8)
    legend("topright", c("Baseline mean", "Scenario mean", "Scenario 95% prediction interval"),
           pch = c(1, 16, NA), lty = c(NA, NA, 1),
           col = c("grey20", "#28658C", "#8BA8BA"), bty = "n", cex = 0.8, bg = "white")
  }
  title(heading, line = 3.2)
  mtext(sprintf("Target %s; CPI %.3f%%, PPI %.3f%% MoM; training cutoff %s",
                config$target_month, config$cpi_mom_pct, config$ppi_mom_pct,
                format(as.Date(cutoff), "%Y-%m")), side = 3, line = 1.3, cex = 0.8)
  unavailable = results$ticker[results$status == "unavailable"]
  missing = if (length(unavailable)) paste(unavailable, collapse = ", ") else "none"
  mtext(sprintf("Modeled %d securities, %.1f%% of security value; unavailable: %s.",
                nrow(d), coverage$modeled_security_coverage_pct, missing),
        side = 1, outer = TRUE, line = 0.5, cex = 0.8)
  mtext(sprintf("Holdings valued %s; cash assumed unchanged; no trades modeled.", as_of),
        side = 1, outer = TRUE, line = 1.6, cex = 0.8)
  mtext(if (kind == "return") "Intervals exclude parameter and lag-selection uncertainty; see table for stationarity flags."
        else "Dollar changes use snapshot values; research estimates with stationarity review flags.",
        side = 1, outer = TRUE, line = 2.7, cex = 0.8)
  invisible(d)
}

save_scenario_plots = function(results, coverage, config, cutoff, as_of) {
  draw = function(kind) plot_scenario_impacts(results, coverage, config, cutoff, as_of, kind)
  grDevices::pdf(file.path(config$output_dir, "scenario_impacts.pdf"), width = 11, height = 9)
  tryCatch({ draw("dollar"); draw("return") }, finally = grDevices::dev.off())
  for (kind in c("return", "dollar")) {
    grDevices::png(file.path(config$output_dir, paste0(kind, "_impact.png")),
                   width = 1650, height = 1350, res = 150)
    tryCatch(draw(kind), finally = grDevices::dev.off())
    if (interactive()) draw(kind)
  }
}

coverage_summary = function(results) {
  modeled = results$status == "modeled"
  total = sum(results$value_usd)
  securities = results$ticker != "CASH"
  security_value = sum(results$value_usd[securities])
  data.frame(total_value_usd = total,
             modeled_value_usd = sum(results$value_usd[modeled]),
             modeled_portfolio_coverage_pct = 100 * sum(results$value_usd[modeled]) / total,
             modeled_security_coverage_pct = if (security_value > 0)
               100 * sum(results$value_usd[modeled]) / security_value else NA_real_,
             cash_value_usd = sum(results$value_usd[!securities]),
             modeled_dollar_difference = if (any(modeled))
               sum(results$dollar_difference[modeled]) else NA_real_)
}

run_cbf_scenarios = function(config = settings, reference = forecast_reference) {
  target = month_date(config$target_month)
  to_log_pct(c(config$cpi_mom_pct, config$ppi_mom_pct))
  holdings = read_holdings(config$holdings_path, config$holdings_as_of)
  config$holdings_path = holdings$path
  config$holdings_as_of = holdings$as_of
  macro = fetch_macro(config$start_date, config$estimation_end_month)
  cutoff = max(macro$month)
  if (target <= cutoff) stop("Target month must follow the common macro cutoff.")
  horizon = length(seq(cutoff, target, by = "month")) - 1L
  if (horizon > 120L) stop("Target is more than 120 months beyond the cutoff.")
  cat("\nMacro cutoff:", format(cutoff, "%Y-%m"), "| Target:", config$target_month,
      "| Holdings valued:", holdings$as_of, "\n")
  cat("Forecast lead:", horizon, "months from the data cutoff.\n")
  cat("Assumptions: CPI", config$cpi_mom_pct, "%, PPI", config$ppi_mom_pct, "% MoM\n")
  cat("Economic-month associations, not release-day reactions or causal effects.\n")
  cat("Returns use Yahoo adjusted prices; dividend treatment is not independently verified.\n")
  config$output_dir = file.path(config$output_dir,
                               paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", config$target_month))
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  results = list(); models = list(); irfs = list(); checks = list()
  p = holdings$positions
  for (j in seq_len(nrow(p))) {
    ticker = if (p$assetType[j] == "cash") "CASH" else p$symbol[j]
    row = data.frame(ticker = ticker, value_usd = p$value_usd[j], sample_months = NA_integer_,
                     selected_lag = NA_integer_, baseline_ppi_pct = NA_real_, baseline_cpi_pct = NA_real_,
                     baseline_return_pct = NA_real_, scenario_return_pct = NA_real_, difference_pp = NA_real_,
                     dollar_difference = NA_real_, lower_95_pct = NA_real_, upper_95_pct = NA_real_,
                     serial_p = NA_real_, status = "unavailable", notes = "")
    if (ticker == "CASH") {
      row$status = "cash assumption"
      row$difference_pp = 0; row$dollar_difference = 0
      row$notes = "Unchanged cash rate; return level not modeled."
    } else {
      cat("Processing", ticker, "...\n")
      row = tryCatch({
        returns = fetch_returns(ticker, config$start_date, cutoff)
        history = aligned_history(macro, returns)
        y = as.matrix(history[, c("PPI", "CPI", "Return")])
        row$sample_months = nrow(y)
        fit = fit_var(y)
        row$selected_lag = fit$p
        row$serial_p = model_checks(fit)
        row$notes = stationarity_notes(y)
        estimates = scenario_forecast(fit, horizon, config$cpi_mom_pct, config$ppi_mom_pct)
        for (name in names(estimates)) row[[name]] = estimates[[name]]
        row$dollar_difference = row$value_usd * row$difference_pp / 100
        row$status = "modeled"
        responses = holding_irfs(fit, config$irf_runs, config$seed)
        bt = backtest(y, config$backtest_months)
        if (nrow(bt)) {
          bt$month = history$month[bt$index]
          bt$ticker = ticker
          checks[[ticker]] = bt
        }
        models[[ticker]] = fit
        irfs[[ticker]] = responses
        row
      }, error = function(e) {
        row$status = "unavailable"
        row$notes = conditionMessage(e)
        for (name in c("baseline_ppi_pct", "baseline_cpi_pct", "baseline_return_pct",
                       "scenario_return_pct", "difference_pp", "dollar_difference",
                       "lower_95_pct", "upper_95_pct")) row[[name]] = NA_real_
        row
      })
    }
    results[[j]] = row
    cat(ticker, ":", row$status, row$notes, "\n")
  }
  table = do.call(rbind, results)
  coverage = coverage_summary(table)
  backtests = if (length(checks)) do.call(rbind, checks) else data.frame()
  accuracy = data.frame()
  if (nrow(backtests)) {
    accuracy = do.call(rbind, lapply(names(checks), function(ticker) {
      b = backtests[backtests$ticker == ticker, ]
      b = b[b$status == "modeled", ]
      data.frame(ticker = ticker, tested_months = nrow(b),
                 var_rmse_pct = if (nrow(b)) sqrt(mean((b$actual_pct - b$var_pct)^2)) else NA_real_,
                 mean_rmse_pct = if (nrow(b)) sqrt(mean((b$actual_pct - b$historical_mean_pct)^2)) else NA_real_,
                 ex_post_conditional_rmse_pct = if (nrow(b))
                   sqrt(mean((b$actual_pct - b$ex_post_conditional_pct)^2)) else NA_real_)
    }))
  }
  write.csv(table, file.path(config$output_dir, "position_scenarios.csv"), row.names = FALSE)
  write.csv(coverage, file.path(config$output_dir, "coverage.csv"), row.names = FALSE)
  write.csv(backtests, file.path(config$output_dir, "backtests.csv"), row.names = FALSE)
  write.csv(accuracy, file.path(config$output_dir, "backtest_accuracy.csv"), row.names = FALSE)
  if (length(models)) {
    lag_selection = do.call(rbind, lapply(names(models), function(ticker) {
      data.frame(ticker = ticker, models[[ticker]]$lag_selection)
    }))
    write.csv(lag_selection, file.path(config$output_dir, "lag_selection.csv"), row.names = FALSE)
  }
  if (!is.null(reference)) write.csv(reference, file.path(config$output_dir, "forecast_reference.csv"), row.names = FALSE)
  if (length(irfs)) {
    grDevices::pdf(file.path(config$output_dir, "impulse_responses.pdf"), width = 11, height = 8)
    tryCatch(for (ticker in names(irfs)) plot_holding_irfs(irfs[[ticker]], ticker, cutoff),
             finally = grDevices::dev.off())
    if (interactive()) for (ticker in names(irfs)) plot_holding_irfs(irfs[[ticker]], ticker, cutoff)
  }
  metadata = list(settings = config, cutoff = as.character(cutoff), holdings_as_of = holdings$as_of,
                  reference = reference,
                  notes = c("Prediction intervals exclude parameter uncertainty.",
                            "Lags minimize BIC among models passing diagnostics; selection uncertainty is excluded.",
                            "Dollar differences use snapshot values, not live valuations.",
                            "Separate VARs have separate macro baselines; no portfolio prediction interval.",
                            "Backtests use revised macro history, not real-time vintages.",
                            "Ex-post conditional backtests use realized inflation, not tradable forecasts.",
                            "IRFs: one standard deviation, 95% pointwise bootstrap bands, months 0-12."))
  saveRDS(list(metadata = metadata, positions = table, coverage = coverage,
               models = models, irfs = irfs, backtests = backtests, accuracy = accuracy),
          file.path(config$output_dir, "research_results.rds"))
  print(table, row.names = FALSE)
  print(coverage, row.names = FALSE)
  if (!is.null(reference)) { cat("\nExisting forecasts (reference only):\n"); print(reference) }
  cat("\n", paste(metadata$notes, collapse = "\n"), "\n", sep = "")
  cat("Results saved to", normalizePath(config$output_dir, winslash = "/"), "\n")
  save_scenario_plots(table, coverage, config, cutoff, holdings$as_of)
  invisible(list(metadata = metadata, positions = table, coverage = coverage,
                 models = models, irfs = irfs, backtests = backtests, accuracy = accuracy))
}

if (isTRUE(getOption("cbf.var.run", TRUE))) cbf_results = run_cbf_scenarios()
