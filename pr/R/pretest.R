library(data.table)

sizeDistortionTest <- function(r, x, lags ) {

  # load tables
  df <- pr::df_gls95

  pretest_table <- pr::pretest

  N <- length(x)


  # Equation (3) of Campbell and Yogo (2005)
  x_diff <- c(NA, diff(x))
  X <- x_diff
  X <- cbind(X, shift(x, n=1, fill=NA, type="lag"))
  if (lags>1){
    for (lag in 1:(lags-1)){
      X <- cbind(X, shift(x_diff, n=lag, fill=NA, type="lag"))

    }
  }
  X <- na.omit(X)
  lm_3 <- lm(X[,1]~ X[,-1])
  smry_3 <- summary(lm_3)
  bic <- BIC(lm_3)

  # Equation (4) of Campbell and Yogo (2005)
  lm_4 <- lm(r[2:N] ~ x[1:(N-1)])
  smry_4 <- summary(lm_4)

  beta_hat <- smry_4$coefficients[2,1]
  SE_beta <- smry_4$coefficients[2,2]


  e <- smry_3$residuals
  u <- smry_4$residuals

  sigma_u_hat_sq <- 1/(length(u)-2)*u%*%u
  sigma_e_hat_sq <- 1/(length(e)-2)*e%*%e
  sigma_ue_hat <- 1/(length(e)-2)*u[(length(u)-length(e)+1):length(u)]%*%e
  delta_hat <- sigma_ue_hat/(sigma_u_hat_sq*sigma_e_hat_sq)^0.5


  #########################DF GLS##########################################

  rho_gls <- 1-7/(length(e))
  y <- x - c(0,rho_gls*x[1:(length(x)-1)])
  x_prime <- c(1,rep(1-rho_gls, length(x)-1))

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
  if (df_gls< (-5)){
    warning("DF-GLS test statistic < -5: Test unreliable.")}

  #get closest tabulated delta value
  table_deltas <- as.numeric(row.names(pretest_table))


  delta_lookup <- table_deltas[which.min(abs(table_deltas - c(delta_hat)) )]


  table_df_gls <- as.numeric(row.names(df))
  df_gls_lookup <- table_df_gls[which.min(abs(table_df_gls - c(df_gls)) )]


  ci <- c(df[as.character(df_gls_lookup),1][[1]], df[as.character(df_gls_lookup),2][[1]])

  ###########test the null hypothesis that the actual size exceeds 7.5%######



  pretest_ci <- c(pretest_table[as.character(delta_lookup),1][[1]],
                  pretest_table[as.character(delta_lookup),2][[1]])


  reject_null <- ( ci[2]<pretest_ci[1] | ci[1]>pretest_ci[2] )

  if (is.na(reject_null)){
    reject_null <- TRUE
    }


  return(reject_null)
}
