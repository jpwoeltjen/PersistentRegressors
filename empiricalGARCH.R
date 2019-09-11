library(rugarch)
require("magrittr")

setwd("~/Desktop/PersistentRegressors")

name <- 'CRSP_M'
dir_name <- paste('data_yogo/',name, '.txt', sep='')
df <- read.csv(dir_name, header=TRUE, sep="\t", dec="." , check.names=TRUE)
r <- df$ret

# standard GARCH gaussian
spec<- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                  mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="norm")
fit <- ugarchfit(spec = spec, data = r)
show(fit)
stargazer::stargazer(fit@fit$matcoef,
                     title = "Parameter Estimates of the GARCH(1, 1)") %>%
  gsub("Std. Error", "Rob. Std. Error", .) %>%
  gsub("t value", "Rob. t value", .) %>%
  # gsub("mu", "$\\\\mu$", .) %>%
  gsub("alpha1", "$\\\\alpha$", .) %>%
  gsub("omega", "$\\\\omega$", .) %>%
  gsub("beta1", "$\\\\beta$", .) %>%
  gsub("gamma1", "$\\\\gamma$", .) %>%
  gsub("shape", "$\\\\theta$", .)  %>%
  writeLines("results/GARCH_output.tex")

# GJR gaussian
spec<- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                  mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="norm")
fit <- ugarchfit(spec = spec, data = r)
show(fit)
stargazer::stargazer(fit@fit$matcoef,
                     title = "Parameter Estimates of the GJR-GARCH(1, 1)") %>%
  gsub("Std. Error", "Rob. Std. Error", .) %>%
  gsub("t value", "Rob. t value", .) %>%
  # gsub("mu", "$\\\\mu$", .) %>%
  gsub("alpha1", "$\\\\alpha$", .) %>%
  gsub("omega", "$\\\\omega$", .) %>%
  gsub("beta1", "$\\\\beta$", .) %>%
  gsub("gamma1", "$\\\\gamma$", .) %>%
  gsub("shape", "$\\\\theta$", .)  %>%
  writeLines("results/GJR-GARCH_norm_output.tex")

# GJR ged
spec<- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                  mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
fit <- ugarchfit(spec = spec, data = r)
show(fit)

stargazer::stargazer(fit@fit$matcoef,
                     title = "Parameter Estimates of the GJR-GARCH(1, 1)") %>%
  gsub("Std. Error", "Rob. Std. Error", .) %>%
  gsub("t value", "Rob. t value", .) %>%
  # gsub("mu", "$\\\\mu$", .) %>%
  gsub("alpha1", "$\\\\alpha$", .) %>%
  gsub("omega", "$\\\\omega$", .) %>%
  gsub("beta1", "$\\\\beta$", .) %>%
  gsub("gamma1", "$\\\\gamma$", .) %>%
  gsub("shape", "$\\\\theta$", .)  %>%
  writeLines("results/GJR-GARCH_ged_output.tex")
