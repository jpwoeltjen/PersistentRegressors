# this file generates mc_output/Cst_GRARCH_table.tex

library(pr)
library(xtable)
library(latex2exp)
library(grDevices)
library(TSA)
library(Matrix)
library(MASS)
library(zoo)
shift <- data.table::shift

Cs <- NULL
rhos <- NULL
deltas <- NULL
obs <- NULL
garch_alphas <- NULL
garch_betas <- NULL
delta_hats <- NULL
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
MC <- 10000
DEBUG <- FALSE
for (N in c(50, 100)){
  N <- N+1 # correct for the loss of 1 observation due to lagging
  for (C in c(-2)){
    print(C)
    for (delta in c(-0.95)){
      for(garch_alpha in c(0.0999)){
        for(garch_beta in c(0.9,0.8,0.5,0.4,0.3,0.2,0.1)){

          sigma_ue <- delta
          Sigma <- matrix(c(sigma_u_sq, sigma_ue, sigma_ue, sigma_e_sq), 2, 2)
          rho <- 1 +C/N

          t_stats <- NULL
          q_bonferroni_two_sided <- NULL
          q_bonferroni_right <- NULL
          infeasible_q <- NULL
          t_right <- NULL
          pretest_rejection <- NULL

          R <- matrix(c(1, delta, delta, 1), 2, 2)
          for (j in 1:MC){
            z <-  mvrnorm(n=N, rep(0, 2), Sigma)
            h <- garch.sim(alpha=c(0.02, garch_alpha), beta=c(garch_beta), n=N)^2
            x <- NULL
            r <- NULL
            x[1] <- 0.001
            for(t in 2:N){
              w = sqrt(h[t])*(z[t,])
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

            t <- s$coefficients[2,3]

            rho_confidence_interval <- rho_ci(x, lags=1, level="0.95")
            t_test_reliable <- as.numeric(sizeDistortionTest(r, x, lags=1))
            ci_beta <- bonferroniQci(r,x, lags=1)
            delta_hat <- ci_beta[5]
            delta_hats <- c(delta_hats,delta_hat)

            q_test_outcome_right_sided <- as.numeric(0<ci_beta[1])
            q_test_outcome_two_sided <- as.numeric((0<ci_beta[1]) | (0>ci_beta[2]))

            if (DEBUG==TRUE){
              print(var(w))
              plot.ts(x)
              plot.ts(r)
              print(s)
              print(t)
              if (ci_beta[7] < (-5)){
                print(ci_beta[7])}
            }


            t_stats[j] <- t
            q_bonferroni_right[j] <- q_test_outcome_right_sided
            t_right[j] <- as.numeric(t > 1.645)
            q_bonferroni_two_sided[j] <- q_test_outcome_two_sided
            pretest_rejection[j] <- t_test_reliable

          }

          Cs <- c(Cs,C)
          rhos <- c(rhos,rho)
          deltas <- c(deltas,delta)
          obs <- c(obs,N-1)
          garch_alphas <- c(garch_alphas,garch_alpha)
          garch_betas <- c(garch_betas,garch_beta)
          bq_rejection_rates <- c(bq_rejection_rates,(mean(q_bonferroni_right)))
          bq_rejection_rates_two <- c(bq_rejection_rates_two, (mean(q_bonferroni_two_sided)))
          t_rejection_rates <- c(t_rejection_rates,(mean(t_right)))
          pretest_rejection_rates <- c(pretest_rejection_rates, (mean(pretest_rejection)))

        }
      }
    }
  }
}
master_df <- data.frame("Obs" = obs,
                        "c" = Cs,
                        "rho" = rhos,
                        "delta" = deltas,
                        "GARCHalpha" = garch_alphas,
                        "GARCHbeta" = garch_betas,
                        "T.test" = t_rejection_rates,
                        "Pretest" = pretest_rejection_rates,
                        "Bonf.Q.test" = bq_rejection_rates,
                        "Bonf.Q.test.two.sided" = bq_rejection_rates_two)

print(xtable(master_df[ , c("Obs","c","rho","delta","GARCHalpha","GARCHbeta","T.test","Bonf.Q.test","Bonf.Q.test.two.sided")],
             digits=c(0,0,0,3,4,4,2,4,4,4), type = "latex"), file = "mc_output/cst_garch_table.tex")
