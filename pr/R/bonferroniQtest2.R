library(data.table)

bonferroniQci <- function(r, x, lags) {


  # load tables
  df <- pr::bonferroniQ

  N <- length(x)

  # Equation (1) of Campbell and Yogo (2005)
  lm_1 <- lm(x[2:N] ~ x[1:(N-1)])
  smry_1 <- summary(lm_1)

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

  rho_hat <- smry_1$coefficients[2,1]
  SE_rho_hat <- smry_1$coefficients[2,2]


  e <- smry_3$residuals
  u <- smry_4$residuals
  v <- smry_1$residuals

  sigma_u_hat_sq <- 1/(length(u)-2)*u%*%u
  sigma_e_hat_sq <- 1/(length(e)-2)*e%*%e
  sigma_ue_hat <- 1/(length(e)-2)*u[(length(u)-length(e)+1):length(u)]%*%e
  sigma_v_hat_sq <- 1/(length(v)-2)*v%*%v
  delta_hat <- sigma_ue_hat/(sigma_u_hat_sq*sigma_e_hat_sq)^0.5
  if (lags>1){
    omega_sq_hat <- sigma_e_hat_sq/(1-sum(lm_3$coefficients[3:length(lm_3$coefficients)]))^2
  }else{
    omega_sq_hat <- sigma_e_hat_sq
  }

  #########################DF GLS##########################################

    rho_gls <- 1-7/(length(e))
    y <- x - c(0,rho_gls*x[1:length(x)-1])
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
      warning("DF-GLS test statistic < -5: confidence intervals unreliable.")}


    #get closest tabulated delta value
    table_deltas <- c(-0.999, -0.975, -0.95, -0.925,
                   -0.9, -0.875, -0.85, -0.825,
                   -0.8, -0.775, -0.75, -0.725,
                   0.7, -0.675, -0.65, -0.625,
                   -0.6, -0.575, -0.55, -0.525,
                   -0.5, -0.475, -0.45, -0.425,
                   -0.4, -0.375, -0.35, -0.325,
                   -0.3, -0.275, -0.25, -0.225,
                   -0.2, -0.175, -0.15, -0.125,
                   -0.1, -0.075, -0.05, -0.025
                   )


  delta_lookup <- table_deltas[which.min(abs(table_deltas - c(delta_hat)) )]


  table_df_gls <- as.numeric(row.names(df))
  df_gls_lookup <- table_df_gls[which.min(abs(table_df_gls - c(df_gls)) )]


  ci <- df[as.character(df_gls_lookup),as.character(delta_lookup)][[1]]

  ci_rho <- c(1+ci[[1]]/length(e), 1+ci[[2]]/length(e))

  r_star_rho_l <- (r[2:length(r)] - c(sigma_ue_hat/(sqrt(sigma_e_hat_sq)*sqrt(omega_sq_hat)))
                   *(x[2:length(x)]-ci_rho[[1]]*x[1:(length(x)-1)]))

  r_star_rho_u <- (r[2:length(r)] - c(sigma_ue_hat/(sqrt(sigma_e_hat_sq)*sqrt(omega_sq_hat)))
                   *(x[2:length(x)]-ci_rho[[2]]*x[1:(length(x)-1)]))

  lm_rstarl <- lm(r_star_rho_l ~ x[1:(length(x)-1)])
  smry_rstarl <- summary(lm_rstarl)
  beta_rho_l <- smry_rstarl$coefficients[2,1]
  se_beta_rho_l <- smry_rstarl$coefficients[2,2]

  lm_rstaru <- lm(r_star_rho_u ~ x[1:(length(x)-1)])
  smry_rstaru <- summary(lm_rstaru)
  beta_rho_u <- smry_rstaru$coefficients[2,1]
  se_beta_rho_u <- smry_rstaru$coefficients[2,2]

  a <- as.numeric((length(e)-2)/2*sigma_ue_hat/(sqrt(sigma_e_hat_sq)*sqrt(omega_sq_hat))
       *(omega_sq_hat/sigma_v_hat_sq-1)*SE_rho_hat^2)

  ci_beta <- c(beta_rho_u+a-1.645*sqrt(1-delta_hat^2)*SE_beta,
               beta_rho_l+a+1.645*sqrt(1-delta_hat^2)*SE_beta)


  return(c(ci_beta, ci_rho, delta_hat, lags, df_gls, sd(e)/sd(u)))

}


