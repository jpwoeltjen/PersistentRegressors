success <- require(urca)

bonferroniQci_ar1 <- function(r,x) {
  print("bonferroniQci_ar1 is deprecated, use bonferroniQci instead!")
  df <- pr::bonferroniQ

  N <- length(x)
  linmod <- lm(r[2:N] ~ x[1:N-1])
  s <- summary(linmod)

  linmod2 <- lm(x[2:N] ~ x[1:N-1])
  s2 <- summary(linmod2)

  beta_hat <- s$coefficients[2,1]
  SE_beta <- s$coefficients[2,2]

  rho_hat <- s2$coefficients[2,1]
  SE_rho <- s2$coefficients[2,2]

  T <- length(x)-1
  u <- s$residuals
  e <- s2$residuals

  sigma_u_hat_sq <- 1/(length(u)-2)*u%*%u
  sigma_e_hat_sq <- 1/(length(e)-2)*e%*%e
  sigma_ue_hat <- 1/(length(e)-2)*u%*%e
  delta_hat <- sigma_ue_hat/(sigma_u_hat_sq*sigma_e_hat_sq)^0.5



  if (success==TRUE){
    ers.x <- ur.ers(x, type="DF-GLS", model="constant", lag.max=0)
    ers.s <- summary(ers.x)
    ers.test_stat <- ers.x@teststat
    df_gls <- ers.test_stat
  }
  else{

    ####DF GLS
    rho_gls <- 1-7/(length(e))
    y <- x -c(0,rho_gls*x[1:length(x)-1])
    x_prime <- c(1,rep(1-rho_gls, length(x)-1))

    # estimate linear model without intercept
    linmod3 <- lm(y~x_prime -1)
    s3 <- summary(linmod3)

    mu_gls <- s3$coefficients[1,1]
    x_bar <-  x - mu_gls
    x_bar_diff <- diff(x_bar)
    x_bar_lag <-  x_bar[1:length(x_bar)-1]

    # estimate linear model without intercept
    linmod4 <- lm( x_bar_diff ~ x_bar_lag -1)
    s4 <- summary(linmod4)
    df_gls <- s4$coefficients["x_bar_lag","t value"]

  }



  table_deltas <- c(-0.999, -0.975, -0.95, -0.925,
                   -0.9, -0.875, -0.85, -0.825,
                   -0.8, -0.775, -0.75, -0.725,
                   0.7, -0.675, -0.65, -0.625,
                   -0.6, -0.575, -0.55, -0.525,
                   -0.5, -0.475, -0.45, -0.425,
                   -0.4, -0.375, -0.35, -0.325,
                   -0.3, -0.275, -0.25, -0.225,
                   -0.2, -0.175, -0.15, -0.125,
                   -0.1, -0.075, -0.05, -0.025)





  delta_lookup <- table_deltas[which.min(abs(table_deltas - c(delta_hat)) )]




  #get closest tabulated df_gls value
  table_df_gls <- as.numeric(row.names(df))
  df_gls_lookup <- table_df_gls[which.min(abs(table_df_gls - c(df_gls)) )]


  ci <- df[as.character(df_gls_lookup),as.character(delta_lookup)][[1]]

  ci_rho <- c(1+ci[[1]]/T, 1+ci[[2]]/T)

  r_star_rho_l <- r[2:length(r)] - c(sigma_ue_hat/sigma_e_hat_sq)*(x[2:length(x)]-ci_rho[[1]]*x[1:(length(x)-1)])
  r_star_rho_u <- r[2:length(r)] - c(sigma_ue_hat/sigma_e_hat_sq)*(x[2:length(x)]-ci_rho[[2]]*x[1:(length(x)-1)])

  linmod5 <- lm(r_star_rho_l ~ x[1:length(x)-1])
  s5 <- summary(linmod5)
  beta_rho_l <- s5$coefficients[2,1]
  se_beta_rho_l <- s5$coefficients[2,2]

  linmod6 <- lm(r_star_rho_u ~ x[1:length(x)-1])
  s6 <- summary(linmod6)
  beta_rho_u <- s6$coefficients[2,1]
  se_beta_rho_u <- s6$coefficients[2,2]

  ci_beta <- c(beta_rho_u-1.645*sqrt(1-delta_hat^2)*SE_beta, beta_rho_l+1.645*sqrt(1-delta_hat^2)*SE_beta)
  return(ci_beta)

}


