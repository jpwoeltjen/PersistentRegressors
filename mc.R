library(pr)
library(xtable)
library(latex2exp)
library(grDevices)

Cs <- NULL
rhos <- NULL
deltas <- NULL
obs <- NULL
t_rejection_rates <- NULL
q_rejection_rates <- NULL
bq_rejection_rates <- NULL
pretest_rejection_rates <- NULL



sigma_u_sq <- 1
sigma_e_sq <- 1
gamma <- 0.0
alpha <- 0.0
beta <- 0.0
beta_0 <- 0
MC <- 10000

for (N in c(50, 10000)){
  N <- N+1 # correct for the loss of 1 observation due to lagging
  for (C in c(0, -2)){
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

        DEBUG=FALSE

        q_test <- function() {
          beta_ue <- sigma_ue/sigma_e_sq
          x_lagged_demeaned <- x[1:length(x)-1] - mean(x[1:length(x)-1])

          return((sum(x_lagged_demeaned*(r[2:length(r)]-beta_0*x[1:length(x)-1]
                                         - beta_ue*(x[2:length(x)]-rho*x[1:length(x)-1])))
                  /((sigma_u_sq*(1-delta^2))^0.5*((sum(x_lagged_demeaned^2))^0.5))
          ))
        }



        for (j in 1:MC){
          w <- mvrnorm(n=N, rep(0, 2), Sigma)


          x <- NULL
          x[1] = w[1,1]
          for(t in 2:N){
            x[t] <- gamma + rho*x[t-1] + w[t,1]
          }

          #AR(2)
          # x <- NULL
          # x[1] = w[1,1]
          # x[2] = w[2,1]
          # for(t in 3:N){
          #   x[t] = gamma + rho*x[t-1] + 0.3*x[t-2] + w[t,1]
          # }

          r <- NULL

          for(t in 2:N){
            r[t] = alpha + beta*x[t-1] + w[t,2]
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
            print(var(w))
            plot.ts(x)
            plot.ts(r)
            print(s)
            print(t)

          }


          t_stats[j] = t
          q_stats[j] = q
          q_bonferroni_right[j] = q_test_outcome_right_sided
          t_right[j] <- as.numeric(t > 1.645)
          q_bonferroni_two_sided[j] = q_test_outcome_two_sided
          infeasible_q[j] = infeasible_q_outcome
          pretest_rejection[j] <- t_test_reliable

        }



        png(paste("mc_output/hist_t_test_", rho,'_' ,delta,'_', N-1, 'Obs','.png', sep=''), width=600, height=600)
        hist(t_stats, density=20, breaks=20, prob=TRUE,
             xlab="t-statistic",
             main=TeX(sprintf("standard normal curve over histogram, $c = %g$, $\\rho = %g$, $\\delta = %g$, $Obs = %g$",
                              C, rho, delta, N-1)))
        curve(dnorm(x, mean=0, sd=1),
              col="darkblue", lwd=2, add=TRUE, yaxt="n")
        dev.off()


        png(paste("mc_output/hist_q_test_", rho,'_', delta, '_', N-1, 'Obs', '.png', sep=''), width=600, height=600)
        hist(q_stats, density=20, breaks=20, prob=TRUE,
             xlab="q-statistic",
             main=TeX(sprintf("standard normal curve over histogram, $c = %g$, $\\rho = %g$, $\\delta = %g$, $Obs = %g$",
                              C, rho, delta, N-1)))
        curve(dnorm(x, mean=0, sd=1),
              col="darkblue", lwd=2, add=TRUE, yaxt="n")
        dev.off()


        Cs <- c(Cs,C)
        rhos <- c(rhos,rho)
        deltas <- c(deltas,delta)
        obs <- c(obs,N-1)

        bq_rejection_rates <- c(bq_rejection_rates,(mean(q_bonferroni_right)))
        # bq_both_rr <- (mean(q_bonferroni_two_sided))
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
                        "Q.test" = q_rejection_rates )

print(xtable(master_df[ , c("Obs","c","rho","delta","T.test","Bonf.Q.test","Q.test")],
             digits=c(0,0,0,3,2,4,4,4), type = "latex"), file = "mc_output/table.tex")
