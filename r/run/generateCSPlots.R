rm(list=ls())

source("../src/models.R")
source("../src/CoxSnellResidualPlot.R")

# Save the first set of plots (Cox-Snell Residual Plot for Cox-PH) as a PNG
png("../figs/CoxSnell_CoxPH_plots.png", width = 1200, height = 800)
par(mfrow = c(1, 3))  # 1 row, 3 columns layout
plotCoxSnellCPH(CoxPH.aic, data.raceCleaned$delta, "AIC")
plotCoxSnellCPH(CoxPH.bic, data.raceCleaned$delta, "BIC")
plotCoxSnellCPH(fullCoxMod, data.raceCleaned$delta, "Full Variable")
dev.off()  # Close the graphics device

# Save the second set of plots (Cox-Snell Residual Plot for Exponential AFT) as a PNG
png("../figs/CoxSnell_ExpAFT_plots.png", width = 1200, height = 800)
par(mfrow = c(1, 3))  # 1 row, 3 columns layout
plotCSExpAFT(expAFT.aic, data.raceCleaned, "AIC")
plotCSExpAFT(expAFT.bic, data.raceCleaned, "BIC")
plotCSExpAFT(full.expAFT, data.raceCleaned, "Full")
dev.off()  # Close the graphics device

# Save the third set of plots (Cox-Snell Residual Plot for Weibull AFT) as a PNG
png("../figs/CoxSnell_WeibullAFT_plots.png", width = 1200, height = 800)
par(mfrow = c(1, 3))  # 1 row, 3 columns layout
plotCSWeibullAFT(weibullAFT.aic, data.raceCleaned, "AIC")
plotCSWeibullAFT(weibullAFT.bic, data.raceCleaned, "BIC")
plotCSWeibullAFT(full.weibullAFT, data.raceCleaned, "Full")
dev.off()  # Close the graphics device

# Save the fourth set of plots (Cox-Snell Residual Plot for Log-Logistic AFT) as a PNG
png("../figs/CoxSnell_LogLogisticAFT_plots.png", width = 1200, height = 800)
par(mfrow = c(1, 3))  # 1 row, 3 columns layout
plotCSLogLogisticAFT(loglogisticAFT.aic, data.raceCleaned, "AIC")
plotCSLogLogisticAFT(loglogisticAFT.bic, data.raceCleaned, "BIC")
plotCSLogLogisticAFT(full.loglogisticAFT, data.raceCleaned, "Full")
dev.off()  # Close the graphics device

# Save the fifth set of plots (Cox-Snell Residual Plot for Log-Normal AFT) as a PNG
png("../figs/CoxSnell_LogNormalAFT_plots.png", width = 1200, height = 800)
par(mfrow = c(1, 3))  # 1 row, 3 columns layout
plotCSLogNormalAFT(lognormalAFT.aic, data.raceCleaned, "AIC")
plotCSLogNormalAFT(lognormalAFT.bic, data.raceCleaned, "BIC")
plotCSLogNormalAFT(full.lognormalAFT, data.raceCleaned, "Full")
dev.off()  # Close the graphics device