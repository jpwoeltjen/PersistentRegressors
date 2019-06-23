\name{rho_ci}
\alias{rho_ci}
\title{Implementing confidence interval for rho according to Campbell and Yogo (2006)}
\usage{
rho_ci(x, lags, level = c("0.95", "0.90", "0.80"), max_lags)
}
\description{
returns ci_rho_lower, ci_rho_upper, ci_c_lower, ci_c_upper, df_gls, opt_lag by BIC
}
\examples{
rho_ci(x, 0, "0.95", max_lags=8)
}
