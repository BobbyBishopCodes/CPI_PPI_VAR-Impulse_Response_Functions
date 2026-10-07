options(cbf.var.run = FALSE)
source("VAR_IRF_CPI_PPI.R")

expect_error = function(expr, pattern) {
  message = tryCatch({ force(expr); "NO ERROR" }, error = function(e) conditionMessage(e))
  stopifnot(grepl(pattern, message))
}
close_to = function(x, y, tolerance = 1e-7) stopifnot(max(abs(x - y)) < tolerance)

close_to(to_simple_pct(to_log_pct(c(-1, 0, 0.3, 2))), c(-1, 0, 0.3, 2))
expect_error(to_log_pct(-100), "greater than")
expect_error(to_log_pct(Inf), "finite")
expect_error(month_date("2026-13"), "month|standard")
expect_error(month_date("October"), "standard")
dates = seq(as.Date("2020-01-01"), by = "month", length.out = 4L)
changes = monthly_changes(dates, c(100, 101, 102, 103), "CPI")
close_to(changes$CPI[1], 100 * log(1.01))
stopifnot(changes$month[1] == dates[2])
expect_error(monthly_changes(dates[-2], c(100, 102, 103), "CPI"), "calendar")
expect_error(monthly_changes(dates, c(100, NA, 102, 103), "CPI"), "invalid")
expect_error(check_months(c(dates[1], dates[1])), "duplicate")
macro = data.frame(month = dates, PPI = 1:4, CPI = 2:5)
returns = data.frame(month = dates[-4], Return = 1:3)
expect_error(aligned_history(macro, returns), "cutoff")

set.seed(744)
n = 600L
A = matrix(c(0.35, 0.12, 0.20, 0, 0.25, -0.15, 0, 0, 0.15), 3, 3)
sigma = matrix(c(0.5, 0.2, 0.3, 0.2, 0.4, -0.1, 0.3, -0.1, 4), 3, 3)
y = matrix(0, n, 3, dimnames = list(NULL, c("PPI", "CPI", "Return")))
innovations = matrix(rnorm(n * 3), n, 3) %*% chol(sigma)
for (j in 2:n) y[j, ] = A %*% y[j - 1L, ] + innovations[j, ]
fit = vars::VAR(y, p = 1L, type = "const")
invisible(model_checks(fit))
selected = fit_var(y)
stopifnot(selected$p >= 1L, selected$p <= 6L)
expect_error(fit_var(y[1:40, ]), "Insufficient history")

d = forecast_distribution(fit, 3L)
S = residual_covariance(fit)
close_to(d$covariance[1:3, 1:3], S)
AA = vars::Acoef(fit)[[1]]
close_to(d$covariance[4:6, 1:3], AA %*% S)
close_to(d$covariance[4:6, 4:6], S + AA %*% S %*% t(AA))
fit2 = vars::VAR(y, p = 2L, type = "const")
d2 = forecast_distribution(fit2, 3L)
coefficients = vars::Acoef(fit2)
phi2 = coefficients[[1]] %*% coefficients[[1]] + coefficients[[2]]
S2 = residual_covariance(fit2)
close_to(d2$covariance[7:9, 1:3], phi2 %*% S2)
conditions = c(0.6, 0.2)
cond = condition_forecast(d, c(7L, 8L), conditions)
close_to(cond$mean[c(7L, 8L)], conditions)
close_to(cond$covariance[c(7L, 8L), ], matrix(0, 2, 9))
stopifnot(min(eigen(cond$covariance, symmetric = TRUE)$values) > -1e-8)
zero = scenario_forecast(fit, 3L, to_simple_pct(d$mean[8L]), to_simple_pct(d$mean[7L]))
close_to(zero[["difference_pp"]], 0)

draws = matrix(rnorm(200000L * 9L), 200000L, 9L) %*% chol(d$covariance)
draws = sweep(draws, 2, d$mean, "+")
gain = t(solve(d$covariance[c(7, 8), c(7, 8)], d$covariance[c(7, 8), ]))
draws = draws + sweep(-draws[, c(7, 8)], 2, conditions, "+") %*% t(gain)
close_to(colMeans(draws), cond$mean, 0.04)
close_to(cov(draws), cond$covariance, 0.05)
r = arithmetic_forecast(cond$mean[9], cond$covariance[9, 9])
close_to(mean(to_simple_pct(draws[, 9])), r[["mean"]], 0.04)
close_to(unname(quantile(to_simple_pct(draws[, 9]), c(0.025, 0.975))),
         unname(r[c("lower", "upper")]), 0.06)

bad = fit
bad$varresult[[1]]$coefficients["PPI.l1"] = 2
expect_error(model_checks(bad), "unstable")
degenerate = vars::VAR(cbind(PPI = y[, 1], CPI = y[, 2], Return = y[, 1]), p = 1, type = "const")
expect_error(residual_covariance(degenerate), "singular")
serial_y = y
for (j in 3:n) serial_y[j, ] = 0.7 * serial_y[j - 2L, ] + innovations[j, ]
expect_error(model_checks(vars::VAR(serial_y, p = 1, type = "const")), "autocorrelation")
serial_fit = fit_var(serial_y)
stopifnot(serial_fit$p >= 2L,
          serial_fit$lag_selection$check[1] != "passed")
invisible(model_checks(serial_fit))
passing = serial_fit$lag_selection[serial_fit$lag_selection$check == "passed", ]
stopifnot(serial_fit$p == passing$lag[which.min(passing$BIC)])

irfs = holding_irfs(fit, runs = 100L, seed = 6L)
irfs_repeat = holding_irfs(fit, runs = 100L, seed = 6L)
stopifnot(identical(irfs, irfs_repeat))
for (impulse in c("PPI", "CPI")) {
  close_to(irfs$PPI_first$irf[[impulse]],
           vars::Psi(fit, nstep = 12L)[3, match(impulse, colnames(y)), ])
}
bt = backtest(y, months = 3L)
stopifnot(nrow(bt) == 3L, all(bt$status == "modeled"))
close_to(bt$historical_mean_pct[1], mean(to_simple_pct(y[1:(n - 3), "Return"])))

rows = data.frame(ticker = c("A", "B", "CASH"), value_usd = c(1000, 2000, 1000),
                  status = c("modeled", "unavailable", "cash assumption"),
                  dollar_difference = c(10, NA, 0))
coverage = coverage_summary(rows)
close_to(coverage$modeled_portfolio_coverage_pct, 25)
close_to(coverage$modeled_security_coverage_pct, 100 / 3)
close_to(coverage$modeled_dollar_difference, 10)
rows$status[1] = "unavailable"
stopifnot(is.na(coverage_summary(rows)$modeled_dollar_difference))

original_fetch_macro = fetch_macro
original_fetch_returns = fetch_returns
synthetic_dates = seq(as.Date("1975-01-01"), by = "month", length.out = n)
fetch_macro = function(start_date, end_month = NULL)
  data.frame(month = synthetic_dates, PPI = y[, "PPI"], CPI = y[, "CPI"])
fetch_returns = function(symbol, start_date, cutoff) {
  if (symbol == "MISSING") stop("Yahoo returned no prices.")
  data.frame(month = synthetic_dates, Return = y[, "Return"])
}
fixture = tempfile(fileext = ".csv")
output = tempfile(pattern = "cbf-var-tests-")
csv = data.frame(Description = c("Example, Stock", "Missing fund", "Bank deposit"),
                 "SYMBOL/CUSIP" = c(" synth ", "MISSING", ""),
                 Quantity = c(5, 10, 1000), "Delayed Price" = c("190", "200", "1.00*"),
                 "Current Value" = c("$1,000.00", "$2,000.00", "$1,000.00"),
                 "Product Type" = c("Stock", "Funds", "Cash & Cash Alternatives"),
                 check.names = FALSE)
write.csv(csv, fixture, row.names = FALSE)
holdings = read_holdings(fixture, "2024-12-31")
stopifnot(identical(holdings$positions$symbol, c("SYNTH", "MISSING", "")),
          identical(holdings$positions$assetType, c("stock", "fund", "cash")),
          holdings$as_of == "2024-12-31")
close_to(holdings$positions$value_usd, c(1000, 2000, 1000))
csv[["Amount Invested (\u2020)"]] = c("900", "1800", "1000")
write.csv(csv, fixture, row.names = FALSE, fileEncoding = "windows-1252")
encoded_holdings = read_holdings(fixture, "2024-12-31")
stopifnot(identical(encoded_holdings$positions, holdings$positions))
write.csv(csv, fixture, row.names = FALSE, fileEncoding = "UTF-8")
bytes = readBin(fixture, "raw", n = file.info(fixture)$size)
connection = file(fixture, "wb")
writeBin(c(as.raw(c(239, 187, 191)), bytes), connection)
close(connection)
stopifnot(identical(read_holdings(fixture, "2024-12-31")$positions, holdings$positions))
invalid = csv
invalid[["SYMBOL/CUSIP"]][2] = "SYNTH"
write.csv(invalid, fixture, row.names = FALSE)
expect_error(read_holdings(fixture, "2024-12-31"), "duplicate")
invalid[["SYMBOL/CUSIP"]][2] = "123456789"
write.csv(invalid, fixture, row.names = FALSE)
expect_error(read_holdings(fixture, "2024-12-31"), "CUSIP")
invalid = csv
invalid[["Current Value"]][1] = "-10"
write.csv(invalid, fixture, row.names = FALSE)
expect_error(read_holdings(fixture, "2024-12-31"), "current values")
write.csv(csv, fixture, row.names = FALSE)
config = settings
config$holdings_path = fixture
config$holdings_as_of = "2024-12-31"
config$start_date = "1975-01-01"
config$target_month = "2025-01"
config$irf_runs = 100L
config$backtest_months = 3L
config$output_dir = output
reference = data.frame(month = "2025-01", CPI = 0.2, PPI = 0.1)
quiet = capture.output(result <- run_cbf_scenarios(config, reference))
stopifnot(nrow(result$positions) == 3L,
          identical(result$positions$status, c("modeled", "unavailable", "cash assumption")),
          is.na(result$positions$scenario_return_pct[2]),
          result$positions$dollar_difference[3] == 0,
          length(result$irfs) == 1L,
          !is.null(result$irfs$SYNTH$PPI_first))
close_to(result$positions$dollar_difference[1], 10 * result$positions$difference_pp[1])
stopifnot(file.exists(file.path(result$metadata$settings$output_dir, "position_scenarios.csv")),
          file.exists(file.path(result$metadata$settings$output_dir, "lag_selection.csv")),
          file.exists(file.path(result$metadata$settings$output_dir, "impulse_responses.pdf")),
          file.exists(file.path(result$metadata$settings$output_dir, "scenario_impacts.pdf")),
          file.exists(file.path(result$metadata$settings$output_dir, "dollar_impact.png")),
          file.exists(file.path(result$metadata$settings$output_dir, "return_impact.png")),
          identical(result$metadata$reference, reference))
config$target_month = "2024-12"
expect_error(run_cbf_scenarios(config), "follow")
fetch_macro = original_fetch_macro
fetch_returns = original_fetch_returns
stopifnot(startsWith(normalizePath(output), normalizePath(tempdir())))
unlink(output, recursive = TRUE)
unlink(fixture)

cat("All VAR/scenario checks passed.\n")
