# univarite GARCH return innovations
library(pr)
library(xtable)
library(latex2exp)
library(grDevices)
library(TSA)
library(Matrix)
library(MASS)
library(zoo)



Cs <- NULL
rhos <- NULL
deltas <- NULL
obs <- NULL
t_rejection_rates <- NULL
q_rejection_rates <- NULL
bq_rejection_rates <- NULL
bq_rejection_rates_two <- NULL
pretest_rejection_rates <- NULL

sigma_u_sq <- 1
sigma_e_sq <- 1
gamma <- 0.0
alpha <- 0.0
beta <- 0.0
beta_0 <- 0
MC <- 1000
DEBUG <- FALSE
for (N in c(100)){
  N <- N+1 # correct for the loss of 1 observation due to lagging
  for (C in c(-2)){
    print(C)
    for (delta in c(-0.95)){

      sigma_ue <- delta
      Sigma <- matrix(c(sigma_u_sq, sigma_ue, sigma_ue, sigma_e_sq), 2, 2)
      rho <- 1 +C/N

      t_stats <- NULL
      q_stats <- NULL
      q_bonferroni_two_sided <- NULL
      q_bonferroni_right <- NULL
      infeasible_q <- NULL
      t_right <- NULL
      pretest_rejection <- NULL

      q_test <- function() {
        beta_ue <- sigma_ue/sigma_e_sq
        x_lagged_demeaned <- x[1:length(x)-1] - mean(x[1:length(x)-1])

        return((sum(x_lagged_demeaned*(r[2:length(r)]-beta_0*x[1:length(x)-1]
                                       - beta_ue*(x[2:length(x)]-rho*x[1:length(x)-1])))
                /((sigma_u_sq*(1-delta^2))^0.5*((sum(x_lagged_demeaned^2))^0.5))
        ))
      }

      garch_alpha <- 0.5
      garch_beta <- 0.4999999
      R <- matrix(c(1, delta, delta, 1), 2, 2)
      for (j in 1:MC){
        z <- mvrnorm(n=N, rep(0, 2), diag(c(1,1)))
        h1 <- garch.sim(alpha=c(1, garch_alpha), beta=c(garch_beta), n=N)



        x <- NULL
        r <- NULL
        x[1] <- 0.001
        for(t in 2:N){
          w = diag(c(1, sqrt(h1_sq[t])))%*%chol(R)%*%(z[t,])
          x[t] <- gamma + rho*x[t-1] + w[1]
          r[t] <- alpha + beta*x[t-1] + w[2]
        }

        linmod <- lm(r[2:N] ~ x[1:N-1])
        s <- summary(linmod)

        linmod2 <- lm(x[2:N] ~ x[1:N-1])
        s2 <- summary(linmod2)

        beta_hat <- s$coefficients[2,1]
        SE_beta <- s$coefficients[2,2]

        rho_hat <- s2$coefficients[2,1]
        SE_rho <- s2$coefficients[2,2]


        std_error_theoretical <- sigma_u_sq^0.5*sum((x[1:N-1]-mean(x[1:N-1]))^2)^-0.5

        t <- s$coefficients[2,3]
        q <- q_test()

        rho_confidence_interval <- rho_ci(x, lags=1, level="0.95")
        t_test_reliable <- as.numeric(sizeDistortionTest(r, x, lags=1))
        ci_beta <- bonferroniQci(r,x, lags=1)


        q_test_outcome_right_sided <- as.numeric(0<ci_beta[1])
        q_test_outcome_two_sided <- as.numeric((0<ci_beta[1]) | (0>ci_beta[2]))
        infeasible_q_outcome <- as.numeric(q > 1.645)

        if (DEBUG==TRUE){
          # print(var(w))
          # plot.ts(x)
          # plot.ts(r)
          # print(s)
          # print(t)
          if (ci_beta[7] < (-5)){
            print(ci_beta[7])}
        }


        t_stats[j] <- t
        q_stats[j] <- q
        q_bonferroni_right[j] <- q_test_outcome_right_sided
        t_right[j] <- as.numeric(t > 1.645)
        q_bonferroni_two_sided[j] <- q_test_outcome_two_sided
        infeasible_q[j] <- infeasible_q_outcome
        pretest_rejection[j] <- t_test_reliable

      }

      Cs <- c(Cs,C)
      rhos <- c(rhos,rho)
      deltas <- c(deltas,delta)
      obs <- c(obs,N-1)

      bq_rejection_rates <- c(bq_rejection_rates,(mean(q_bonferroni_right)))
      bq_rejection_rates_two <- c(bq_rejection_rates_two, (mean(q_bonferroni_two_sided)))
      q_rejection_rates <- c(q_rejection_rates,(mean(infeasible_q)))
      t_rejection_rates <- c(t_rejection_rates,(mean(t_right)))
      pretest_rejection_rates <- c(pretest_rejection_rates, (mean(pretest_rejection)))

    }
  }
}

master_df <- data.frame("Obs" = obs,
                        "c" = Cs,
                        "rho" = rhos,
                        "delta" = deltas,
                        "T.test" = t_rejection_rates,
                        "Pretest" = pretest_rejection_rates,
                        "Bonf.Q.test" = bq_rejection_rates,
                        "Bonf.Q.test.two.sided" = bq_rejection_rates_two)

