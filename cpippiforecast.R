#CPI and PPI Forecasting - Riley's Rewritten Version
#Target: September 2026 | CPI release 10/14/2026 | PPI release 10/15/2026

# ============================================================
# LIBRARIES
# ============================================================
library(quantmod)
library(xts)
library(zoo)
library(tseries)
library(TSA)
library(vars)
library(forecast)

# ============================================================
# PART 1: DATA PREPARATION AND VISUALIZATION
# ============================================================

# --- Data Import ---
start_date = "2010-01-01"
end_date = "2026-08-31"

getSymbols("CPIAUCSL", src = "FRED")
getSymbols("PPIFIS", src = "FRED")

# --- Trim to the sample window (FRED pulls the full series) ---
cpi_level = CPIAUCSL[paste0(start_date, "/", end_date)]
ppi_level = PPIFIS[paste0(start_date, "/", end_date)]

# --- Rename columns cleanly ---
colnames(cpi_level) = "CPI"
colnames(ppi_level) = "PPI"

# --- Align dates across both series (inner join) ---
all_data <- na.omit(merge(cpi_level, ppi_level))

# --- Convert to ts objects (monthly, starting Jan 2010) ---
start_yr <- as.numeric(format(index(all_data)[1], "%Y"))
start_mo <- as.numeric(format(index(all_data)[1], "%m"))

CPI = ts(as.numeric(all_data$CPI), start = c(start_yr, start_mo), frequency = 12)
PPI = ts(as.numeric(all_data$PPI), start = c(start_yr, start_mo), frequency = 12)

# --- Individual plots ---
plot(CPI, col='red')
plot(PPI, col='blue')

# --- Multiple plots
plot(CPI, col='red', xlab = 'Time', ylab = 'Index')
lines(PPI, col='blue')

# CPI is the headline consumer price index, seasonally adjusted
# PPI is the final demand producer price index, seasonally adjusted

# --- Computations
mean(CPI)
median(CPI)
var(CPI)

mean(PPI)
median(PPI)
var(PPI)

prices = cbind(CPI, PPI)
cor(prices)

# --- Variable of Interest (VoI)
plot(CPI, col='red')

# --- Models of VoI
m1 = lm(CPI~time(CPI))
summary(m1)

m2 = lm(CPI~time(CPI)+I(time(CPI)^2))
summary(m2)

# --- Actual and Fitted Values
f1=ts(m1$fitted.values,start = start(CPI),
      frequency = frequency(CPI))
plot(CPI, col='red')
lines(f1, col='blue')

f2=ts(m2$fitted.values,start = start(CPI),
      frequency = frequency(CPI))
plot(CPI, col='red')
lines(f2, col='darkgreen')

# --- Diagnostic Check
res1=rstudent(m1)
plot(res1,x=as.vector(time(CPI)),type='o')

res2=rstudent(m2)
plot(res2,x=as.vector(time(CPI)),type='o')

hist(res1)
hist(res2)

qqnorm(res1)
qqline(res1)
qqnorm(res2)
qqline(res2)

shapiro.test(res1)
shapiro.test(res2)

Box.test(res1, lag = 12, type = 'Ljung')
Box.test(res2, lag = 12, type = 'Ljung')

acf(res1)
acf(res2)

# --- Comparing Models
summary(m1)
summary(m2)

AIC(m1)
AIC(m2)

# --- Stationarity and Transformations
# The fund forecasts the month-over-month percent change, which is what BLS
# headlines, so the index gets differenced into exact percent changes.
plot(CPI, col='red')

adf.test(CPI)
kpss.test(CPI, null='Level')
kpss.test(CPI, null='Trend')

CPIg = 100 * (CPI / stats::lag(CPI, -1) - 1)
PPIg = 100 * (PPI / stats::lag(PPI, -1) - 1)

adf.test(CPIg)
kpss.test(CPIg, null='Level')
kpss.test(CPIg, null='Trend')

adf.test(PPIg)
kpss.test(PPIg, null='Level')
kpss.test(PPIg, null='Trend')

plot(CPIg, col='red', ylab = 'CPI MoM (%)')
plot(PPIg, col='blue', ylab = 'PPI MoM (%)')

# ============================================================
# PART 2: CPI MODEL SELECTION
# ============================================================

# --- ARIMA Model Specification
acf(CPIg)
pacf(CPIg)
eacf(CPIg)

# --- ARIMA Estimation and Diagnostics
cmodel1 = arima(CPIg, order=c(0,0,0))
cmodel1

cmodel2 = arima(CPIg, order=c(1,0,0))
cmodel2

cmodel3 = arima(CPIg, order=c(1,0,1))
cmodel3

cmodel4 = arima(CPIg, order=c(1,0,2))
cmodel4

AIC(cmodel1,cmodel2,cmodel3,cmodel4)

# --- Diagnostics
cres1=residuals(cmodel1)
cres2=residuals(cmodel2)
cres3=residuals(cmodel3)
cres4=residuals(cmodel4)

plot(cres4)
hist(cres4)
acf(cres4)

qqnorm(cres4)
qqline(cres4)
shapiro.test(cres4)

LB.test(cmodel1)
LB.test(cmodel2)
LB.test(cmodel3)
LB.test(cmodel4)

# ============================================================
# PART 3: PPI MODEL SELECTION
# ============================================================

# --- ARIMA Model Specification
acf(PPIg)
pacf(PPIg)
eacf(PPIg)

# --- ARIMA Estimation and Diagnostics
pmodel1 = arima(PPIg, order=c(0,0,0))
pmodel1

pmodel2 = arima(PPIg, order=c(1,0,0))
pmodel2

pmodel3 = arima(PPIg, order=c(1,0,1))
pmodel3


AIC(pmodel1,pmodel2,pmodel3)

# --- Diagnostics
pres1=residuals(pmodel1)
pres2=residuals(pmodel2)
pres3=residuals(pmodel3)

plot(pres3)
hist(pres3)
acf(pres3)

qqnorm(pres3)
qqline(pres3)
shapiro.test(pres3)

LB.test(pmodel1)
LB.test(pmodel2)
LB.test(pmodel3)

# ============================================================
# PART 4: CPI FORECASTING
# ============================================================
# ARIMA(1,0,2) was the preferred specification in Part 2, so it is the only
# one carried into the out of sample check.

# --- insample OOsample
series = window(CPIg, end = c(2025, 8))
future = window(CPIg, start = c(2025, 9))

acf(series)
pacf(series)
eacf(series)

h = length(future)
n = length(series)

model4 = arima(series, order = c(1, 0, 2))
model4
pred4 = predict(model4, n.ahead = length(future))
pred4$pred
plot(series, col='red', xlim = c(2023, 2027), main = "CPI ARIMA(1,0,2)")
lines(pred4$pred, col='blue', lwd = 2)
lines(future, col='gray30', lwd = 2)

mse4 = mean((as.numeric(future) - as.numeric(pred4$pred))^2)
mse4

AIC(model4)

# ============================================================
# PART 5: PPI FORECASTING
# ============================================================
# ARIMA(1,0,1) was the preferred specification in Part 3.

# --- insample OOsample
pseries = window(PPIg, end = c(2025, 8))
pfuture = window(PPIg, start = c(2025, 9))

acf(pseries)
pacf(pseries)
eacf(pseries)

pmod3 = arima(pseries, order = c(1, 0, 1))
pmod3
ppred3 = predict(pmod3, n.ahead = length(pfuture))
ppred3$pred
plot(pseries, col='blue', xlim = c(2023, 2027), main = "PPI ARIMA(1,0,1)")
lines(ppred3$pred, col='royalblue3', lwd = 2)
lines(pfuture, col='gray30', lwd = 2)

pmse3 = mean((as.numeric(pfuture) - as.numeric(ppred3$pred))^2)
pmse3

AIC(pmod3)

# ============================================================
# PART 6: SEPTEMBER 2026 FORECAST
# ============================================================
# Refit the selected order on the full sample and step one month ahead.

cpi_final = arima(CPIg, order = c(1, 0, 2))
cpi_final
cpi_fc = predict(cpi_final, n.ahead = 1)
cpi_point = as.numeric(cpi_fc$pred)
cpi_lower = cpi_point - 1.2816 * as.numeric(cpi_fc$se)
cpi_upper = cpi_point + 1.2816 * as.numeric(cpi_fc$se)

ppi_final = arima(PPIg, order = c(1, 0, 1))
ppi_final
ppi_fc = predict(ppi_final, n.ahead = 1)
ppi_point = as.numeric(ppi_fc$pred)
ppi_lower = ppi_point - 1.2816 * as.numeric(ppi_fc$se)
ppi_upper = ppi_point + 1.2816 * as.numeric(ppi_fc$se)

# --- CPI forecast plot
plot(CPIg, col='red', xlim = c(2024, 2027), ylab = 'Change from previous month (%)',
     main = "CPI Forecast - September 2026")
abline(h = 0, lty = 3, col='gray70')
lines(cpi_fc$pred, col='royalblue3', type = 'o', pch = 19, lwd = 2)
arrows(2026 + 8/12, cpi_lower, 2026 + 8/12, cpi_upper,
       angle = 90, code = 3, length = 0.07, col='royalblue3', lwd = 2)
legend('topleft', c('Actual','Forecast','80% range'), col = c('red','royalblue3','royalblue3'),
       lty = c(1,1,1), pch = c(NA,19,NA), bty = 'n', cex = 0.8)

# --- PPI forecast plot
plot(PPIg, col='blue', xlim = c(2024, 2027), ylab = 'Change from previous month (%)',
     main = "PPI Forecast - September 2026")
abline(h = 0, lty = 3, col='gray70')
lines(ppi_fc$pred, col='royalblue3', type = 'o', pch = 19, lwd = 2)
arrows(2026 + 8/12, ppi_lower, 2026 + 8/12, ppi_upper,
       angle = 90, code = 3, length = 0.07, col='royalblue3', lwd = 2)
legend('topleft', c('Actual','Forecast','80% range'), col = c('blue','royalblue3','royalblue3'),
       lty = c(1,1,1), pch = c(NA,19,NA), bty = 'n', cex = 0.8)

# ============================================================
# PART 7: MULTIVARIATE ANALYSIS
# ============================================================
# Does producer price inflation lead consumer price inflation?

inflation = na.omit(cbind(CPIg, PPIg))
colnames(inflation) = c("CPI", "PPI")

vs = VARselect(inflation, lag.max = 6, type = "const")
vs

optimal_lag = vs$selection["AIC(n)"]
optimal_lag

var_const = VAR(inflation, p = optimal_lag, type = "const")
var_none  = VAR(inflation, p = optimal_lag, type = "none")
summary(var_const)

causality(var_const, cause = "PPI")
causality(var_const, cause = "CPI")

irf1 = irf(var_const, impulse = "PPI", response = "CPI",
           n.ahead = 12, boot = TRUE)
plot(irf1)

irf2 = irf(var_const, impulse = "CPI", response = "PPI",
           n.ahead = 12, boot = TRUE)
plot(irf2)

# ============================================================
# RESULTS
# ============================================================
cat("CPI Sep 2026:", round(cpi_point, 2), "% | 80% range",
    round(cpi_lower, 2), "to", round(cpi_upper, 2), "| ARIMA(1,0,2) | release 2026-10-14\n")
cat("PPI Sep 2026:", round(ppi_point, 2), "% | 80% range",
    round(ppi_lower, 2), "to", round(ppi_upper, 2), "| ARIMA(1,0,1) | release 2026-10-15\n")

