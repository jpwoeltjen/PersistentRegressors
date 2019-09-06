library(pr)
library(zoo)
library(data.table)
library(ggplot2)
library(xtable)


setwd("~/Desktop/PersistentRegressors")

name <- 'SP_A'
dir_name <- paste('data_yogo/',name, '.txt', sep='')
df <- read.csv(dir_name, header=TRUE, sep="\t", dec="." , check.names=TRUE)

max_year <- 2003
min_year <- 1880

name_list <- NULL
var_list <- NULL
delta_list <- NULL
df_gls_list <- NULL
lags_list <- NULL
t_stat_list <- NULL
pretest_list <- NULL
beta_list <- NULL
rho_ci_list <- NULL
beta_ci_scaled_list <- NULL
obs_list <- NULL




############## MA(10) of earnings / Index price ###############################

df2<- na.omit(df[,c('time','ret', 'lep')])
df2 <- df2[df2$time<=max_year, ]
df2 <- df2[df2$time>=min_year, ]

x <- df2$lep
r <- df2$ret

lm1 <- lm(r[-1] ~ x[-length(x)])

smry1 <- summary(lm1)
beta = smry1$coefficients[2,1]
t_stat = smry1$coefficients[2,3]

rho_ci_outcome <- rho_ci(x, lags=0, level="0.95", max_lags=6)
opt_lags <- rho_ci_outcome[6]
rho_ci_bonferroni <- rho_ci_outcome[1:2]
c_ci_bonferroni <- rho_ci_outcome[3:4]

t_test_reliable <- sizeDistortionTest(r, x, lags=opt_lags)

bonferroniQ_outcome <- bonferroniQci(r, x, lags=opt_lags)

beta_ci <- bonferroniQ_outcome[1:2]
beta_ci_scaled <- round(beta_ci*bonferroniQ_outcome[8],3)
rho_ci_ref_bonferroni <- round(bonferroniQ_outcome[3:4],3)
delta_hat <- round(bonferroniQ_outcome[5],3)
df_gls <- round(bonferroniQ_outcome[7],3)

name_list <- c(name_list, name)
var_list <- c(var_list, colnames(df2)[3])
delta_list <- c(delta_list, delta_hat)
df_gls_list <- c(df_gls_list, df_gls)
lags_list <- c(lags_list, opt_lags)
t_stat_list <- c(t_stat_list, t_stat)
pretest_list <- c(pretest_list, t_test_reliable)
beta_list <- c(beta_list, beta)
rho_ci_list <-c(rho_ci_list, paste('[',rho_ci_ref_bonferroni[1],',',rho_ci_ref_bonferroni[2],']', sep=''))
beta_ci_scaled_list <-c(beta_ci_scaled_list, paste('[',beta_ci_scaled[1],',',beta_ci_scaled[2],']', sep=''))
obs_list <- c(obs_list, length(x))


############## 1 yyyy dividends / Index price ###############################

df2<- na.omit(df[,c('time','ret', 'ldp')])
df2 <- df2[df2$time<=max_year,]
df2 <- df2[df2$time>=min_year,]


x <- df2$ldp
r <- df2$ret

lm1 <- lm(r[-1] ~ x[-length(x)])
smry1 <- summary(lm1)
beta = smry1$coefficients[2,1]
t_stat = smry1$coefficients[2,3]

rho_ci_outcome <- rho_ci(x, lags=0, level="0.95", max_lags=6)
opt_lags <- rho_ci_outcome[6]
rho_ci_bonferroni <- rho_ci_outcome[1:2]
c_ci_bonferroni <- rho_ci_outcome[3:4]

t_test_reliable <- sizeDistortionTest(r, x, lags=opt_lags)

bonferroniQ_outcome <- bonferroniQci(r, x, lags=opt_lags)

beta_ci <- bonferroniQ_outcome[1:2]
beta_ci_scaled <- round(beta_ci*bonferroniQ_outcome[8],3)
rho_ci_ref_bonferroni <- round(bonferroniQ_outcome[3:4],3)
delta_hat <- round(bonferroniQ_outcome[5],3)
df_gls <- round(bonferroniQ_outcome[7],3)

name_list <- c(name_list, name)
var_list <- c(var_list, colnames(df2)[3])
delta_list <- c(delta_list, delta_hat)
df_gls_list <- c(df_gls_list, df_gls)
lags_list <- c(lags_list, opt_lags)
t_stat_list <- c(t_stat_list, t_stat)
pretest_list <- c(pretest_list, t_test_reliable)
beta_list <- c(beta_list, beta)
rho_ci_list <-c(rho_ci_list, paste('[',rho_ci_ref_bonferroni[1],',',rho_ci_ref_bonferroni[2],']', sep=''))
beta_ci_scaled_list <-c(beta_ci_scaled_list, paste('[',beta_ci_scaled[1],',',beta_ci_scaled[2],']', sep=''))
obs_list <- c(obs_list, length(x))

master_df <- data.frame(  name_list,
                          var_list,
                          obs_list,
                          delta_list,
                          rho_ci_list,
                          df_gls_list,
                          lags_list,
                          t_stat_list,
                          as.numeric(pretest_list),
                          beta_list,
                          beta_ci_scaled_list)

# file_name <- paste('results/yogo_empirical_results_', Sys.time(), '.tex', sep='')
# print(xtable(master_df, digits=c(0, 0,0,0,2,0,3,0,2,0,3,0), type = "latex"), file = file_name, include.rownames=FALSE)

