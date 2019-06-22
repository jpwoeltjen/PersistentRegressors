library(data.table)

rho_ci <- function(r, x, lags, level = c("0.95", "0.90", "0.80")) {

  level <- match.arg(level)

  # load tables
  switch(level,
         "0.95" = df <- pr::df_gls95,
         "0.90" = df <- pr::df_gls90,
         "0.80"  = df <- pr::df_gls80)


  N <- length(x)

  rho_gls <- 1-7/(N-1)
  y <- x - c(0,rho_gls*x[1:(N-1)])
  x_prime <- c(1,rep(1-rho_gls, (N-1)))

  # estimate linear model without intercept
  lm_mugls <- lm(y~x_prime -1)
  smry_mugls <- summary(lm_mugls)

  mu_gls <- smry_mugls$coefficients[1,1]
  x_bar <-  x - mu_gls
  x_bar_diff <- c(NA, diff(x_bar))

  X <- x_bar_diff
  X <- cbind(X, shift(x_bar, n=1, fill=NA, type="lag"))
  if (lags>1){
    for (lag in 1:(lags-1)){
      X <- cbind(X, shift(x_bar_diff, n=lag, fill=NA, type="lag"))

    }
  }
  X <- na.omit(X)

  # estimate linear model without intercept
  lm_dfgls <- lm( X[,1]~ X[,-1]  -1)
  smry_dfgls <- summary(lm_dfgls)
  df_gls <- smry_dfgls$coefficients[1,"t value"]

  table_df_gls <- as.numeric(row.names(df))
  df_gls_lookup <- table_df_gls[which.min(abs(table_df_gls - c(df_gls)) )]


  ci <- c(df[as.character(df_gls_lookup),1][[1]], df[as.character(df_gls_lookup),2][[1]])
  ci_rho <- c(1+ci[[1]]/(N-1), 1+ci[[2]]/(N-1))

  return(c(ci_rho, ci))
}
