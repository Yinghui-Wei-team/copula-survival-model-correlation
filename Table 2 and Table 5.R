#######################################################################################
#
#   Filename    :	      Table 2 and Table 5.R    												  
#                                                                                     
#   Project     :       Article "Estimating the correlation between semi-competing risk survival endpoints"                                                             
#   Authors     :       L Sorrell, Y Wei, M Wojtys and P Rowe                                                               
#   Date        :       01/06/2021
#																				  
#   R Version   :       R-3.6.1                                                              
#
#   Required R packages :  copula, mvtnorm, mstate
#
########################################################################################
library(copula) 
library(mvtnorm)
library(mstate) #ACS data
source("Functions.R")   #Likelihood functions

data(aidssi2) #Load ACS data set

X <- aidssi2$si.time-aidssi2$entry.time #time to non-terminal event (SI switch)
d1 <- aidssi2$si.stat #indicator for SI switch
Y <- aidssi2$death.time-aidssi2$entry.time #time to terminal event (death from AIDS)
d2 <- aidssi2$death.stat #indicator for death from AIDS
df <- data.frame(X,Y,d1,d2)

########################
## Recreating Table 2 ## 
########################

table2 <- addmargins(table(d2,d1))

########################
## Recreating Table 5 ## 
########################

#prepare table as data frame:
lambda1 <- rep(NA, 4) #lambda_1
lambda1_lwci <- rep(NA, 4) #lower confidence interval for lambda_1
lambda1_upci <- rep(NA,4) #upper confidence interval for lambda_1
lambda2 <- rep(NA, 4) #lambda_2
lambda2_lwci <- rep(NA, 4) #lower confidence interval for lambda_2
lambda2_upci <- rep(NA,4) #upper confidence interval for lambda_2
rho <- rep(NA, 4) #rho
rho_lwci <- rep(NA, 4) #lower confidence interval for rho
rho_upci <- rep(NA,4) #upper confidence interval for rho
AIC <- rep(NA,4) #AIC
copula <- c("Normal", "Clayton", "Frank", "Gumbel")
table_5 <- data.frame(copula, lambda1, lambda1_lwci, lambda1_upci, lambda2, lambda2_lwci, lambda2_upci, rho, rho_lwci, rho_upci, AIC)


##Clayton & Exponential copula model##
clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.1),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

clayton_fisher_info <- solve(-clayton_optim$hessian) 
clayton_se <- sqrt(diag(clayton_fisher_info)) 

#95% CI for hazard rates and association parameter
clayton_upper_ci_l1 <- clayton_optim$par[1]+1.96*clayton_se[1]
clayton_lower_ci_l1 <- clayton_optim$par[1]-1.96*clayton_se[1]
clayton_upper_ci_l2 <- clayton_optim$par[2]+1.96*clayton_se[2]
clayton_lower_ci_l2 <- clayton_optim$par[2]-1.96*clayton_se[2]
clayton_upper_ci_theta <- clayton_optim$par[3]+1.96*clayton_se[3]
clayton_lower_ci_theta <- clayton_optim$par[3]-1.96*clayton_se[3]

#95% CI for association parameter to rho
clayton_upper_ci_theta_copula <- claytonCopula(clayton_upper_ci_theta)
clayton_upper_ci_rho <- rho(clayton_upper_ci_theta_copula)
clayton_lower_ci_theta_copula <- claytonCopula(clayton_lower_ci_theta)
clayton_lower_ci_rho <- rho(clayton_lower_ci_theta_copula)

clayton_l1 <- clayton_optim$par[1]
clayton_l2 <- clayton_optim$par[2]
clayton_theta <- clayton_optim$par[3]
clayton_copula <- claytonCopula(clayton_theta)
clayton_rho <- rho(clayton_copula)

#AIC
clayton_estimated_parameters <- c(clayton_l1, clayton_l2, clayton_theta)
clayton_loglik_estimated_parameters <- clayton_loglik(clayton_estimated_parameters,X,Y,d1,d2)
k <- length(clayton_estimated_parameters)
clayton_aic <- -2*clayton_loglik_estimated_parameters+2*k

#add to table 5
table_5$lambda1[2] <- clayton_l1
table_5$lambda1_lwci[2] <- clayton_lower_ci_l1
table_5$lambda1_upci[2] <- clayton_upper_ci_l1
table_5$lambda2[2] <- clayton_l2
table_5$lambda2_lwci[2] <- clayton_lower_ci_l2
table_5$lambda2_upci[2] <- clayton_upper_ci_l2
table_5$rho[2] <- clayton_rho
table_5$rho_lwci[2] <- clayton_lower_ci_rho
table_5$rho_upci[2] <- clayton_upper_ci_rho
table_5$AIC[2] <- clayton_aic


##Frank & Exponential copula analysis##
frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.1),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

frank_fisher_info <- solve(-frank_optim$hessian) 
frank_se <- sqrt(diag(frank_fisher_info)) 

#95% CI for hazard rates and association parameter
frank_upper_ci_l1 <- frank_optim$par[1]+1.96*frank_se[1]
frank_lower_ci_l1 <- frank_optim$par[1]-1.96*frank_se[1]
frank_upper_ci_l2 <- frank_optim$par[2]+1.96*frank_se[2]
frank_lower_ci_l2 <- frank_optim$par[2]-1.96*frank_se[2]
frank_upper_ci_theta <- frank_optim$par[3]+1.96*frank_se[3]
frank_lower_ci_theta <- frank_optim$par[3]-1.96*frank_se[3]

#95% CI for association parameter to rho
frank_upper_ci_theta_copula <- frankCopula(frank_upper_ci_theta) #convert CI to rho
frank_upper_ci_rho <- rho(frank_upper_ci_theta_copula)
frank_lower_ci_theta_copula <- frankCopula(frank_lower_ci_theta)
frank_lower_ci_rho <- rho(frank_lower_ci_theta_copula)

frank_l1 <- frank_optim$par[1]
frank_l2 <- frank_optim$par[2]
frank_theta <- frank_optim$par[3]
frank_copula <- frankCopula(frank_theta)
frank_rho <- rho(frank_copula)

#AIC
frank_estimated_parameters <- c(frank_l1, frank_l2, frank_theta)
frank_loglik_estimated_parameters <- frank_loglik(frank_estimated_parameters,X,Y,d1,d2)
k <- length(frank_estimated_parameters)
frank_aic <- -2*frank_loglik_estimated_parameters+2*k

#add to table 5
table_5$lambda1[3] <- frank_l1
table_5$lambda1_lwci[3] <- frank_lower_ci_l1
table_5$lambda1_upci[3] <- frank_upper_ci_l1
table_5$lambda2[3] <- frank_l2
table_5$lambda2_lwci[3] <- frank_lower_ci_l2
table_5$lambda2_upci[3] <- frank_upper_ci_l2
table_5$rho[3] <- frank_rho
table_5$rho_lwci[3] <- frank_lower_ci_rho
table_5$rho_upci[3] <- frank_upper_ci_rho
table_5$AIC[3] <- frank_aic


##Gumbel & Exponential copula analysis##
gumbel_optim <- optim(c(0.1,0.1,0.4), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.01,1.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

gumbel_fisher_info <- solve(-gumbel_optim$hessian) 
gumbel_se <- sqrt(diag(gumbel_fisher_info)) 

#95% CI for hazard rates and association parameter
gumbel_upper_ci_l1 <- gumbel_optim$par[1]+1.96*gumbel_se[1]
gumbel_lower_ci_l1 <- gumbel_optim$par[1]-1.96*gumbel_se[1]
gumbel_upper_ci_l2 <- gumbel_optim$par[2]+1.96*gumbel_se[2]
gumbel_lower_ci_l2 <- gumbel_optim$par[2]-1.96*gumbel_se[2]
gumbel_upper_ci_theta <- gumbel_optim$par[3]+1.96*gumbel_se[3]
gumbel_lower_ci_theta <- gumbel_optim$par[3]-1.96*gumbel_se[3]

#95% CI for association parameter to rho
gumbel_upper_ci_theta_copula <- gumbelCopula(gumbel_upper_ci_theta) #convert CI to rho
gumbel_upper_ci_rho <- rho(gumbel_upper_ci_theta_copula)
gumbel_lower_ci_theta_copula <- gumbelCopula(gumbel_lower_ci_theta)
gumbel_lower_ci_rho <- rho(gumbel_lower_ci_theta_copula)

gumbel_l1 <- gumbel_optim$par[1]
gumbel_l2 <- gumbel_optim$par[2]
gumbel_theta <- gumbel_optim$par[3]
gumbel_copula <- gumbelCopula(gumbel_theta)
gumbel_rho <- rho(gumbel_copula)

#AIC
gumbel_estimated_parameters <- c(gumbel_l1, gumbel_l2, gumbel_theta)
gumbel_loglik_estimated_parameters <- gumbel_loglik(gumbel_estimated_parameters,X,Y,d1,d2)
k <- length(gumbel_estimated_parameters)
gumbel_aic <- -2*gumbel_loglik_estimated_parameters+2*k

#add to table 5
table_5$lambda1[4] <- gumbel_l1
table_5$lambda1_lwci[4] <- gumbel_lower_ci_l1
table_5$lambda1_upci[4] <- gumbel_upper_ci_l1
table_5$lambda2[4] <- gumbel_l2
table_5$lambda2_lwci[4] <- gumbel_lower_ci_l2
table_5$lambda2_upci[4] <- gumbel_upper_ci_l2
table_5$rho[4] <- gumbel_rho
table_5$rho_lwci[4] <- gumbel_lower_ci_rho
table_5$rho_upci[4] <- gumbel_upper_ci_rho
table_5$AIC[4] <- gumbel_aic


##Normal & Exponential copula analysis##
normal_optim <- optim(c(0.1,0.1,0.4), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

normal_fisher_info <- solve(-normal_optim$hessian) 
normal_se <- sqrt(diag(normal_fisher_info)) 

#95% CI for hazard rates and association parameter
normal_upper_ci_l1 <- normal_optim$par[1]+1.96*normal_se[1]
normal_lower_ci_l1 <- normal_optim$par[1]-1.96*normal_se[1]
normal_upper_ci_l2 <- normal_optim$par[2]+1.96*normal_se[2]
normal_lower_ci_l2 <- normal_optim$par[2]-1.96*normal_se[2]
normal_upper_ci_theta <- normal_optim$par[3]+1.96*normal_se[3]
normal_lower_ci_theta <- normal_optim$par[3]-1.96*normal_se[3]

#95% CI for association parameter to rho
normal_upper_ci_theta_copula <- normalCopula(normal_upper_ci_theta) #convert CI to rho
normal_upper_ci_rho <- rho(normal_upper_ci_theta_copula)
normal_lower_ci_theta_copula <- normalCopula(normal_lower_ci_theta)
normal_lower_ci_rho <- rho(normal_lower_ci_theta_copula)

normal_l1 <- normal_optim$par[1]
normal_l2 <- normal_optim$par[2]
normal_theta <- normal_optim$par[3]
normal_copula <- normalCopula(normal_theta)
normal_rho <- rho(normal_copula)

#AIC
normal_estimated_parameters <- c(normal_l1, normal_l2, normal_theta)
normal_loglik_estimated_parameters <- normal_loglik(normal_estimated_parameters,X, Y, d1, d2)
k <- length(normal_estimated_parameters)
normal_aic <- -2*normal_loglik_estimated_parameters+2*k

#add to table 5
table_5$lambda1[1] <- normal_l1
table_5$lambda1_lwci[1] <- normal_lower_ci_l1
table_5$lambda1_upci[1] <- normal_upper_ci_l1
table_5$lambda2[1] <- normal_l2
table_5$lambda2_lwci[1] <- normal_lower_ci_l2
table_5$lambda2_upci[1] <- normal_upper_ci_l2
table_5$rho[1] <- normal_rho
table_5$rho_lwci[1] <- normal_lower_ci_rho
table_5$rho_upci[1] <- normal_upper_ci_rho
table_5$AIC[1] <- normal_aic

#########################
## Recreating Table S3 ## 
#########################

#prepare table as data frame:
alpha1 <- rep(NA, 4) #alpha_1
alpha1_lwci <- rep(NA, 4) #lower confidence interval for alpha_1
alpha1_upci <- rep(NA,4) #upper confidence interval for alpha_1
beta1 <- rep(NA, 4) #beta_1
beta1_lwci <- rep(NA, 4) #lower confidence interval for beta_1
beta1_upci <- rep(NA,4) #upper confidence interval for beta_1
alpha2 <- rep(NA, 4) #alpha_2
alpha2_lwci <- rep(NA, 4) #lower confidence interval for alpha_2
alpha2_upci <- rep(NA,4) #upper confidence interval for alpha_2
beta2 <- rep(NA, 4) #beta_2
beta2_lwci <- rep(NA, 4) #lower confidence interval for beta_2
beta2_upci <- rep(NA,4) #upper confidence interval for beta_2
rho <- rep(NA, 4) #rho
rho_lwci <- rep(NA, 4) #lower confidence interval for rho
rho_upci <- rep(NA,4) #upper confidence interval for rho
AIC <- rep(NA,4) #AIC
table_S3 <- data.frame(alpha1, alpha1_lwci, alpha1_upci, beta1, beta1_lwci, beta1_upci, 
                       alpha2, alpha2_lwci, alpha2_upci, beta2, beta2_lwci, beta2_upci, 
                       rho, rho_lwci, rho_upci, AIC)


##Clayton & Weibull copula analysis##
clayton_weibull_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_weibull_loglik, method="L-BFGS-B",
                               lower=c(0.01,0.01,0.01,0.0001,0.01),upper=c(1.5,1,2,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

clayton_weibull_fisher_info <- solve(-clayton_weibull_optim$hessian) 
clayton_weibull_se <- sqrt(diag(clayton_weibull_fisher_info)) 

#95% CI for hazard rates and association parameter
clayton_weibull_upper_ci_a1 <- clayton_weibull_optim$par[1]+1.96*clayton_weibull_se[1]
clayton_weibull_lower_ci_a1 <- clayton_weibull_optim$par[1]-1.96*clayton_weibull_se[1]
clayton_weibull_upper_ci_b1 <- clayton_weibull_optim$par[2]+1.96*clayton_weibull_se[2]
clayton_weibull_lower_ci_b1 <- clayton_weibull_optim$par[2]-1.96*clayton_weibull_se[2]
clayton_weibull_upper_ci_a2 <- clayton_weibull_optim$par[3]+1.96*clayton_weibull_se[3]
clayton_weibull_lower_ci_a2 <- clayton_weibull_optim$par[3]-1.96*clayton_weibull_se[3]
clayton_weibull_upper_ci_b2 <- clayton_weibull_optim$par[4]+1.96*clayton_weibull_se[4]
clayton_weibull_lower_ci_b2 <- clayton_weibull_optim$par[4]-1.96*clayton_weibull_se[4]
clayton_weibull_upper_ci_theta <- clayton_weibull_optim$par[5]+1.96*clayton_weibull_se[5]
clayton_weibull_lower_ci_theta <- clayton_weibull_optim$par[5]-1.96*clayton_weibull_se[5]

#95% CI for association parameter to rho
clayton_weibull_upper_ci_theta_copula <- claytonCopula(clayton_weibull_upper_ci_theta) #convert CI to rho
clayton_weibull_upper_ci_rho <- rho(clayton_weibull_upper_ci_theta_copula)
clayton_weibull_lower_ci_theta_copula <- claytonCopula(clayton_weibull_lower_ci_theta)
clayton_weibull_lower_ci_rho <- rho(clayton_weibull_lower_ci_theta_copula)

clayton_weibull_a1 <- clayton_weibull_optim$par[1]
clayton_weibull_b1 <- clayton_weibull_optim$par[2]
clayton_weibull_a2 <- clayton_weibull_optim$par[3]
clayton_weibull_b2 <- clayton_weibull_optim$par[4]
clayton_weibull_theta <- clayton_weibull_optim$par[5]
clayton_weibull_copula <- claytonCopula(clayton_weibull_theta)
clayton_weibull_rho <- rho(clayton_weibull_copula)

#AIC
clayton_weibull_estimated_parameters <- c(clayton_weibull_a1, clayton_weibull_b1, clayton_weibull_a2, clayton_weibull_b2, clayton_weibull_theta)
clayton_weibull_loglik_estimated_parameters <- clayton_weibull_loglik(clayton_weibull_estimated_parameters,X,Y,d1,d2)
k <- length(clayton_weibull_estimated_parameters)
clayton_weibull_aic <- -2*clayton_weibull_loglik_estimated_parameters+2*k

table_S3$alpha1[2] <- clayton_weibull_a1
table_S3$alpha1_lwci[2] <- clayton_weibull_lower_ci_a1
table_S3$alpha1_upci[2] <- clayton_weibull_upper_ci_a1
table_S3$beta1[2] <- clayton_weibull_b1
table_S3$beta1_lwci[2] <- clayton_weibull_lower_ci_b1
table_S3$beta1_upci[2] <- clayton_weibull_upper_ci_b1
table_S3$alpha2[2] <- clayton_weibull_a2
table_S3$alpha2_lwci[2] <- clayton_weibull_lower_ci_a2
table_S3$alpha2_upci[2] <- clayton_weibull_upper_ci_a2
table_S3$beta2[2] <- clayton_weibull_b2
table_S3$beta2_lwci[2] <- clayton_weibull_lower_ci_b2
table_S3$beta2_upci[2] <- clayton_weibull_upper_ci_b2
table_S3$rho[2] <- clayton_weibull_rho
table_S3$rho_lwci[2] <- clayton_weibull_lower_ci_rho
table_S3$rho_upci[2] <- clayton_weibull_upper_ci_rho
table_S3$AIC[2] <- clayton_weibull_aic


##Frank & Weibull copula analysis##
frank_weibull_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik, method="L-BFGS-B",
                             lower=c(0.01,0.01,0.01,0.0001,0.01),upper=c(1.5,1,2,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

frank_weibull_fisher_info <- solve(-frank_weibull_optim$hessian) 
frank_weibull_se <- sqrt(diag(frank_weibull_fisher_info)) 

#95% CI for hazard rates and association parameter
frank_weibull_upper_ci_a1 <- frank_weibull_optim$par[1]+1.96*frank_weibull_se[1]
frank_weibull_lower_ci_a1 <- frank_weibull_optim$par[1]-1.96*frank_weibull_se[1]
frank_weibull_upper_ci_b1 <- frank_weibull_optim$par[2]+1.96*frank_weibull_se[2]
frank_weibull_lower_ci_b1 <- frank_weibull_optim$par[2]-1.96*frank_weibull_se[2]
frank_weibull_upper_ci_a2 <- frank_weibull_optim$par[3]+1.96*frank_weibull_se[3]
frank_weibull_lower_ci_a2 <- frank_weibull_optim$par[3]-1.96*frank_weibull_se[3]
frank_weibull_upper_ci_b2 <- frank_weibull_optim$par[4]+1.96*frank_weibull_se[4]
frank_weibull_lower_ci_b2 <- frank_weibull_optim$par[4]-1.96*frank_weibull_se[4]
frank_weibull_upper_ci_theta <- frank_weibull_optim$par[5]+1.96*frank_weibull_se[5]
frank_weibull_lower_ci_theta <- frank_weibull_optim$par[5]-1.96*frank_weibull_se[5]

#95% CI for association parameter to rho
frank_weibull_upper_ci_theta_copula <- frankCopula(frank_weibull_upper_ci_theta) #convert CI to rho
frank_weibull_upper_ci_rho <- rho(frank_weibull_upper_ci_theta_copula)
frank_weibull_lower_ci_theta_copula <- frankCopula(frank_weibull_lower_ci_theta)
frank_weibull_lower_ci_rho <- rho(frank_weibull_lower_ci_theta_copula)

frank_weibull_a1 <- frank_weibull_optim$par[1]
frank_weibull_b1 <- frank_weibull_optim$par[2]
frank_weibull_a2 <- frank_weibull_optim$par[3]
frank_weibull_b2 <- frank_weibull_optim$par[4]
frank_weibull_theta <- frank_weibull_optim$par[5]
frank_weibull_copula <- frankCopula(frank_weibull_theta)
frank_weibull_rho <- rho(frank_weibull_copula)

#AIC
frank_weibull_estimated_parameters <- c(frank_weibull_a1, frank_weibull_b1, frank_weibull_a2, frank_weibull_b2, frank_weibull_theta)
frank_weibull_loglik_estimated_parameters <- frank_weibull_loglik(frank_weibull_estimated_parameters,X,Y,d1,d2)
k <- length(frank_weibull_estimated_parameters)
frank_weibull_aic <- -2*frank_weibull_loglik_estimated_parameters+2*k

table_S3$alpha1[3] <- frank_weibull_a1
table_S3$alpha1_lwci[3] <- frank_weibull_lower_ci_a1
table_S3$alpha1_upci[3] <- frank_weibull_upper_ci_a1
table_S3$beta1[3] <- frank_weibull_b1
table_S3$beta1_lwci[3] <- frank_weibull_lower_ci_b1
table_S3$beta1_upci[3] <- frank_weibull_upper_ci_b1
table_S3$alpha2[3] <- frank_weibull_a2
table_S3$alpha2_lwci[3] <- frank_weibull_lower_ci_a2
table_S3$alpha2_upci[3] <- frank_weibull_upper_ci_a2
table_S3$beta2[3] <- frank_weibull_b2
table_S3$beta2_lwci[3] <- frank_weibull_lower_ci_b2
table_S3$beta2_upci[3] <- frank_weibull_upper_ci_b2
table_S3$rho[3] <- frank_weibull_rho
table_S3$rho_lwci[3] <- frank_weibull_lower_ci_rho
table_S3$rho_upci[3] <- frank_weibull_upper_ci_rho
table_S3$AIC[3] <- frank_weibull_aic

##Gumbel & Weibull copula analysis##
gumbel_weibull_optim <- optim(c(0.2,0.2,0.2,0.2,2), gumbel_weibull_loglik, method="L-BFGS-B",
                              lower=c(0.01,0.01,0.01,0.00001,1.01),upper=c(1,1,2,2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

gumbel_weibull_fisher_info <- solve(-gumbel_weibull_optim$hessian) 
gumbel_weibull_se <- sqrt(diag(gumbel_weibull_fisher_info)) 

#95% CI for hazard rates and association parameter
gumbel_weibull_upper_ci_a1 <- gumbel_weibull_optim$par[1]+1.96*gumbel_weibull_se[1]
gumbel_weibull_lower_ci_a1 <- gumbel_weibull_optim$par[1]-1.96*gumbel_weibull_se[1]
gumbel_weibull_upper_ci_b1 <- gumbel_weibull_optim$par[2]+1.96*gumbel_weibull_se[2]
gumbel_weibull_lower_ci_b1 <- gumbel_weibull_optim$par[2]-1.96*gumbel_weibull_se[2]
gumbel_weibull_upper_ci_a2 <- gumbel_weibull_optim$par[3]+1.96*gumbel_weibull_se[3]
gumbel_weibull_lower_ci_a2 <- gumbel_weibull_optim$par[3]-1.96*gumbel_weibull_se[3]
gumbel_weibull_upper_ci_b2 <- gumbel_weibull_optim$par[4]+1.96*gumbel_weibull_se[4]
gumbel_weibull_lower_ci_b2 <- gumbel_weibull_optim$par[4]-1.96*gumbel_weibull_se[4]
gumbel_weibull_upper_ci_theta <- gumbel_weibull_optim$par[5]+1.96*gumbel_weibull_se[5]
gumbel_weibull_lower_ci_theta <- gumbel_weibull_optim$par[5]-1.96*gumbel_weibull_se[5]

#95% CI for association parameter to rho
gumbel_weibull_upper_ci_theta_copula <- gumbelCopula(gumbel_weibull_upper_ci_theta) #convert CI to rho
gumbel_weibull_upper_ci_rho <- rho(gumbel_weibull_upper_ci_theta_copula)
gumbel_weibull_lower_ci_theta_copula <- gumbelCopula(gumbel_weibull_lower_ci_theta)
gumbel_weibull_lower_ci_rho <- rho(gumbel_weibull_lower_ci_theta_copula)

gumbel_weibull_a1 <- gumbel_weibull_optim$par[1]
gumbel_weibull_b1 <- gumbel_weibull_optim$par[2]
gumbel_weibull_a2 <- gumbel_weibull_optim$par[3]
gumbel_weibull_b2 <- gumbel_weibull_optim$par[4]
gumbel_weibull_theta <- gumbel_weibull_optim$par[5]
gumbel_weibull_copula <- gumbelCopula(gumbel_weibull_theta)
gumbel_weibull_rho <- rho(gumbel_weibull_copula)

#AIC
gumbel_weibull_estimated_parameters <- c(gumbel_weibull_a1, gumbel_weibull_b1, gumbel_weibull_a2, gumbel_weibull_b2, gumbel_weibull_theta)
gumbel_weibull_loglik_estimated_parameters <- gumbel_weibull_loglik(gumbel_weibull_estimated_parameters,X,Y,d1,d2)
k <- length(gumbel_weibull_estimated_parameters)
gumbel_weibull_aic <- -2*gumbel_weibull_loglik_estimated_parameters+2*k

table_S3$alpha1[4] <- gumbel_weibull_a1
table_S3$alpha1_lwci[4] <- gumbel_weibull_lower_ci_a1
table_S3$alpha1_upci[4] <- gumbel_weibull_upper_ci_a1
table_S3$beta1[4] <- gumbel_weibull_b1
table_S3$beta1_lwci[4] <- gumbel_weibull_lower_ci_b1
table_S3$beta1_upci[4] <- gumbel_weibull_upper_ci_b1
table_S3$alpha2[4] <- gumbel_weibull_a2
table_S3$alpha2_lwci[4] <- gumbel_weibull_lower_ci_a2
table_S3$alpha2_upci[4] <- gumbel_weibull_upper_ci_a2
table_S3$beta2[4] <- gumbel_weibull_b2
table_S3$beta2_lwci[4] <- gumbel_weibull_lower_ci_b2
table_S3$beta2_upci[4] <- gumbel_weibull_upper_ci_b2
table_S3$rho[4] <- gumbel_weibull_rho
table_S3$rho_lwci[4] <- gumbel_weibull_lower_ci_rho
table_S3$rho_upci[4] <- gumbel_weibull_upper_ci_rho
table_S3$AIC[4] <- gumbel_weibull_aic

##Normal & Weibull copula analysis##
normal_weibull_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik, method="L-BFGS-B",
                              lower=c(0.01,0.01,0.01,0.00001,0.1),upper=c(1.5,0.2,2,0.1,0.9), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

normal_weibull_fisher_info <- solve(-normal_weibull_optim$hessian) 
normal_weibull_se <- sqrt(diag(normal_weibull_fisher_info)) 

#95% CI for hazard rates and association parameter
normal_weibull_upper_ci_a1 <- normal_weibull_optim$par[1]+1.96*normal_weibull_se[1]
normal_weibull_lower_ci_a1 <- normal_weibull_optim$par[1]-1.96*normal_weibull_se[1]
normal_weibull_upper_ci_b1 <- normal_weibull_optim$par[2]+1.96*normal_weibull_se[2]
normal_weibull_lower_ci_b1 <- normal_weibull_optim$par[2]-1.96*normal_weibull_se[2]
normal_weibull_upper_ci_a2 <- normal_weibull_optim$par[3]+1.96*normal_weibull_se[3]
normal_weibull_lower_ci_a2 <- normal_weibull_optim$par[3]-1.96*normal_weibull_se[3]
normal_weibull_upper_ci_b2 <- normal_weibull_optim$par[4]+1.96*normal_weibull_se[4]
normal_weibull_lower_ci_b2 <- normal_weibull_optim$par[4]-1.96*normal_weibull_se[4]
normal_weibull_upper_ci_theta <- normal_weibull_optim$par[5]+1.96*normal_weibull_se[5]
normal_weibull_lower_ci_theta <- normal_weibull_optim$par[5]-1.96*normal_weibull_se[5]

#95% CI for association parameter to rho
normal_weibull_upper_ci_theta_copula <- normalCopula(normal_weibull_upper_ci_theta) #convert CI to rho
normal_weibull_upper_ci_rho <- rho(normal_weibull_upper_ci_theta_copula)
normal_weibull_lower_ci_theta_copula <- normalCopula(normal_weibull_lower_ci_theta)
normal_weibull_lower_ci_rho <- rho(normal_weibull_lower_ci_theta_copula)

normal_weibull_a1 <- normal_weibull_optim$par[1]
normal_weibull_b1 <- normal_weibull_optim$par[2]
normal_weibull_a2 <- normal_weibull_optim$par[3]
normal_weibull_b2 <- normal_weibull_optim$par[4]
normal_weibull_theta <- normal_weibull_optim$par[5]
normal_weibull_copula <- normalCopula(normal_weibull_theta)
normal_weibull_rho <- rho(normal_weibull_copula)

#AIC
normal_weibull_estimated_parameters <- c(normal_weibull_a1, normal_weibull_b1, normal_weibull_a2, normal_weibull_b2, normal_weibull_theta)
normal_weibull_loglik_estimated_parameters <- normal_weibull_loglik(normal_weibull_estimated_parameters,X,Y,d1,d2)
k <- length(normal_weibull_estimated_parameters)
normal_weibull_aic <- -2*normal_weibull_loglik_estimated_parameters+2*k

table_S3$alpha1[1] <- normal_weibull_a1
table_S3$alpha1_lwci[1] <- normal_weibull_lower_ci_a1
table_S3$alpha1_upci[1] <- normal_weibull_upper_ci_a1
table_S3$beta1[1] <- normal_weibull_b1
table_S3$beta1_lwci[1] <- normal_weibull_lower_ci_b1
table_S3$beta1_upci[1] <- normal_weibull_upper_ci_b1
table_S3$alpha2[1] <- normal_weibull_a2
table_S3$alpha2_lwci[1] <- normal_weibull_lower_ci_a2
table_S3$alpha2_upci[1] <- normal_weibull_upper_ci_a2
table_S3$beta2[1] <- normal_weibull_b2
table_S3$beta2_lwci[1] <- normal_weibull_lower_ci_b2
table_S3$beta2_upci[1] <- normal_weibull_upper_ci_b2
table_S3$rho[1] <- normal_weibull_rho
table_S3$rho_lwci[1] <- normal_weibull_lower_ci_rho
table_S3$rho_upci[1] <- normal_weibull_upper_ci_rho
table_S3$AIC[1] <- normal_weibull_aic


#########################
## Recreating Table S4 ## 
#########################

#prepare table as data frame:
gamma1 <- rep(NA, 4) #gamma_1
gamma1_lwci <- rep(NA, 4) #lower confidence interval for gamma_1
gamma1_upci <- rep(NA,4) #upper confidence interval for gamma_1
lambda1 <- rep(NA, 4) #lambda_1
lambda1_lwci <- rep(NA, 4) #lower confidence interval for lambda_1
lambda1_upci <- rep(NA,4) #upper confidence interval for lambda_1
gamma2 <- rep(NA, 4) #gamma_2
gamma2_lwci <- rep(NA, 4) #lower confidence interval for gamma_2
gamma2_upci <- rep(NA,4) #upper confidence interval for gamma_2
lambda2 <- rep(NA, 4) #lambda_2
lambda2_lwci <- rep(NA, 4) #lower confidence interval for lambda_2
lambda2_upci <- rep(NA,4) #upper confidence interval for lambda_2
rho <- rep(NA, 4) #rho
rho_lwci <- rep(NA, 4) #lower confidence interval for rho
rho_upci <- rep(NA,4) #upper confidence interval for rho
AIC <- rep(NA,4) #AIC
table_S4 <- data.frame(gamma1, gamma1_lwci, gamma1_upci, lambda1, lambda1_lwci, lambda1_upci, 
                     gamma2, gamma2_lwci, gamma2_upci, lambda2, lambda2_lwci, lambda2_upci, 
                     rho, rho_lwci, rho_upci, AIC)


##Clayton & Gompertz copula analysis##
clayton_gompertz_optim <- optim(c(0.1,0.1,0.1,0.1,2), clayton_gompertz_loglik, method="L-BFGS-B",
                                lower=c(-0.1,0.01,-0.1,0.01,0.01),upper=c(0.1,0.1,0.6,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

clayton_gompertz_fisher_info <- solve(-clayton_gompertz_optim$hessian) 
clayton_gompertz_se <- sqrt(diag(clayton_gompertz_fisher_info)) 

#95% CI for hazard rates and association parameter
clayton_gompertz_upper_ci_g1 <- clayton_gompertz_optim$par[1]+1.96*clayton_gompertz_se[1]
clayton_gompertz_lower_ci_g1 <- clayton_gompertz_optim$par[1]-1.96*clayton_gompertz_se[1]
clayton_gompertz_upper_ci_l1 <- clayton_gompertz_optim$par[2]+1.96*clayton_gompertz_se[2]
clayton_gompertz_lower_ci_l1 <- clayton_gompertz_optim$par[2]-1.96*clayton_gompertz_se[2]
clayton_gompertz_upper_ci_g2 <- clayton_gompertz_optim$par[3]+1.96*clayton_gompertz_se[3]
clayton_gompertz_lower_ci_g2 <- clayton_gompertz_optim$par[3]-1.96*clayton_gompertz_se[3]
clayton_gompertz_upper_ci_l2 <- clayton_gompertz_optim$par[4]+1.96*clayton_gompertz_se[4]
clayton_gompertz_lower_ci_l2 <- clayton_gompertz_optim$par[4]-1.96*clayton_gompertz_se[4]
clayton_gompertz_upper_ci_theta <- clayton_gompertz_optim$par[5]+1.96*clayton_gompertz_se[5]
clayton_gompertz_lower_ci_theta <- clayton_gompertz_optim$par[5]-1.96*clayton_gompertz_se[5]

#95% CI for association parameter to rho
clayton_gompertz_upper_ci_theta_copula <- claytonCopula(clayton_gompertz_upper_ci_theta) #convert CI to rho
clayton_gompertz_upper_ci_rho <- rho(clayton_gompertz_upper_ci_theta_copula)
clayton_gompertz_lower_ci_theta_copula <- claytonCopula(clayton_gompertz_lower_ci_theta)
clayton_gompertz_lower_ci_rho <- rho(clayton_gompertz_lower_ci_theta_copula)

clayton_gompertz_g1 <- clayton_gompertz_optim$par[1]
clayton_gompertz_l1 <- clayton_gompertz_optim$par[2]
clayton_gompertz_g2 <- clayton_gompertz_optim$par[3]
clayton_gompertz_l2 <- clayton_gompertz_optim$par[4]
clayton_gompertz_theta <- clayton_gompertz_optim$par[5]
clayton_gompertz_copula <- claytonCopula(clayton_gompertz_theta)
clayton_gompertz_rho <- rho(clayton_gompertz_copula)

#AIC
clayton_gompertz_estimated_parameters <- c(clayton_gompertz_g1, clayton_gompertz_l1, clayton_gompertz_g2, clayton_gompertz_l2, clayton_gompertz_theta)
clayton_gompertz_loglik_estimated_parameters <- clayton_gompertz_loglik(clayton_gompertz_estimated_parameters,X,Y,d1,d2)
k <- length(clayton_gompertz_estimated_parameters)
clayton_gompertz_aic <- -2*clayton_gompertz_loglik_estimated_parameters+2*k

table_S4$gamma1[2] <- clayton_gompertz_g1
table_S4$gamma1_lwci[2] <- clayton_gompertz_lower_ci_g1
table_S4$gamma1_upci[2] <- clayton_gompertz_upper_ci_g1
table_S4$lambda1[2] <- clayton_gompertz_l1
table_S4$lambda1_lwci[2] <- clayton_gompertz_lower_ci_l1
table_S4$lambda1_upci[2] <- clayton_gompertz_upper_ci_l1
table_S4$gamma2[2] <-clayton_gompertz_g2
table_S4$gamma2_lwci[2] <- clayton_gompertz_lower_ci_g2
table_S4$gamma2_upci[2] <- clayton_gompertz_upper_ci_g2
table_S4$lambda2[2] <- clayton_gompertz_l2
table_S4$lambda2_lwci[2] <- clayton_gompertz_lower_ci_l2
table_S4$lambda2_upci[2] <- clayton_gompertz_upper_ci_l2
table_S4$rho[2] <- clayton_gompertz_rho
table_S4$rho_lwci[2] <- clayton_gompertz_lower_ci_rho
table_S4$rho_upci[2] <- clayton_gompertz_upper_ci_rho
table_S4$AIC[2] <- clayton_gompertz_aic


##Frank & Gompertz copula analysis##
frank_gompertz_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik, method="L-BFGS-B",
                              lower=c(-0.1,0.01,-0.1,0.01,0.01),upper=c(0.1,0.1,0.6,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

frank_gompertz_fisher_info <- solve(-frank_gompertz_optim$hessian) 
frank_gompertz_se <- sqrt(diag(frank_gompertz_fisher_info)) 

#95% CI for hazard rates and association parameter
frank_gompertz_upper_ci_g1 <- frank_gompertz_optim$par[1]+1.96*frank_gompertz_se[1]
frank_gompertz_lower_ci_g1 <- frank_gompertz_optim$par[1]-1.96*frank_gompertz_se[1]
frank_gompertz_upper_ci_l1 <- frank_gompertz_optim$par[2]+1.96*frank_gompertz_se[2]
frank_gompertz_lower_ci_l1 <- frank_gompertz_optim$par[2]-1.96*frank_gompertz_se[2]
frank_gompertz_upper_ci_g2 <- frank_gompertz_optim$par[3]+1.96*frank_gompertz_se[3]
frank_gompertz_lower_ci_g2 <- frank_gompertz_optim$par[3]-1.96*frank_gompertz_se[3]
frank_gompertz_upper_ci_l2 <- frank_gompertz_optim$par[4]+1.96*frank_gompertz_se[4]
frank_gompertz_lower_ci_l2 <- frank_gompertz_optim$par[4]-1.96*frank_gompertz_se[4]
frank_gompertz_upper_ci_theta <- frank_gompertz_optim$par[5]+1.96*frank_gompertz_se[5]
frank_gompertz_lower_ci_theta <- frank_gompertz_optim$par[5]-1.96*frank_gompertz_se[5]

#95% CI for association parameter to rho
frank_gompertz_upper_ci_theta_copula <- frankCopula(frank_gompertz_upper_ci_theta) #convert CI to rho
frank_gompertz_upper_ci_rho <- rho(frank_gompertz_upper_ci_theta_copula)
frank_gompertz_lower_ci_theta_copula <- frankCopula(frank_gompertz_lower_ci_theta)
frank_gompertz_lower_ci_rho <- rho(frank_gompertz_lower_ci_theta_copula)

frank_gompertz_g1 <- frank_gompertz_optim$par[1]
frank_gompertz_l1 <- frank_gompertz_optim$par[2]
frank_gompertz_g2 <- frank_gompertz_optim$par[3]
frank_gompertz_l2 <- frank_gompertz_optim$par[4]
frank_gompertz_theta <- frank_gompertz_optim$par[5]
frank_gompertz_copula <- frankCopula(frank_gompertz_theta)
frank_gompertz_rho <- rho(frank_gompertz_copula)

#AIC
frank_gompertz_estimated_parameters <- c(frank_gompertz_g1, frank_gompertz_l1, frank_gompertz_g2, frank_gompertz_l2, frank_gompertz_theta)
frank_gompertz_loglik_estimated_parameters <- frank_gompertz_loglik(frank_gompertz_estimated_parameters,X,Y,d1,d2)
k <- length(frank_gompertz_estimated_parameters)
frank_gompertz_aic <- -2*frank_gompertz_loglik_estimated_parameters+2*k

table_S4$gamma1[3] <- frank_gompertz_g1
table_S4$gamma1_lwci[3] <- frank_gompertz_lower_ci_g1
table_S4$gamma1_upci[3] <- frank_gompertz_upper_ci_g1
table_S4$lambda1[3] <- frank_gompertz_l1
table_S4$lambda1_lwci[3] <- frank_gompertz_lower_ci_l1
table_S4$lambda1_upci[3] <- frank_gompertz_upper_ci_l1
table_S4$gamma2[3] <- frank_gompertz_g2
table_S4$gamma2_lwci[3] <- frank_gompertz_lower_ci_g2
table_S4$gamma2_upci[3] <- frank_gompertz_upper_ci_g2
table_S4$lambda2[3] <- frank_gompertz_l2
table_S4$lambda2_lwci[3] <- frank_gompertz_lower_ci_l2
table_S4$lambda2_upci[3] <- frank_gompertz_upper_ci_l2
table_S4$rho[3] <- frank_gompertz_rho
table_S4$rho_lwci[3] <- frank_gompertz_lower_ci_rho
table_S4$rho_upci[3] <- frank_gompertz_upper_ci_rho
table_S4$AIC[3] <- frank_gompertz_aic


##Gumbel & Gompertz copula analysis##
gumbel_gompertz_optim <- optim(c(0.1,0.1,0.1,0.1,2), gumbel_gompertz_loglik, method="L-BFGS-B",
                               lower=c(-0.1,0.01,-0.1,0.01,1.01),upper=c(0.1,0.1,0.6,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

gumbel_gompertz_fisher_info <- solve(-gumbel_gompertz_optim$hessian) 
gumbel_gompertz_se <- sqrt(diag(gumbel_gompertz_fisher_info)) 

#95% CI for hazard rates and association parameter
gumbel_gompertz_upper_ci_g1 <- gumbel_gompertz_optim$par[1]+1.96*gumbel_gompertz_se[1]
gumbel_gompertz_lower_ci_g1 <- gumbel_gompertz_optim$par[1]-1.96*gumbel_gompertz_se[1]
gumbel_gompertz_upper_ci_l1 <- gumbel_gompertz_optim$par[2]+1.96*gumbel_gompertz_se[2]
gumbel_gompertz_lower_ci_l1 <- gumbel_gompertz_optim$par[2]-1.96*gumbel_gompertz_se[2]
gumbel_gompertz_upper_ci_g2 <- gumbel_gompertz_optim$par[3]+1.96*gumbel_gompertz_se[3]
gumbel_gompertz_lower_ci_g2 <- gumbel_gompertz_optim$par[3]-1.96*gumbel_gompertz_se[3]
gumbel_gompertz_upper_ci_l2 <- gumbel_gompertz_optim$par[4]+1.96*gumbel_gompertz_se[4]
gumbel_gompertz_lower_ci_l2 <- gumbel_gompertz_optim$par[4]-1.96*gumbel_gompertz_se[4]
gumbel_gompertz_upper_ci_theta <- gumbel_gompertz_optim$par[5]+1.96*gumbel_gompertz_se[5]
gumbel_gompertz_lower_ci_theta <- gumbel_gompertz_optim$par[5]-1.96*gumbel_gompertz_se[5]

#95% CI for association parameter to rho
gumbel_gompertz_upper_ci_theta_copula <- gumbelCopula(gumbel_gompertz_upper_ci_theta) #convert CI to rho
gumbel_gompertz_upper_ci_rho <- rho(gumbel_gompertz_upper_ci_theta_copula)
gumbel_gompertz_lower_ci_theta_copula <- gumbelCopula(gumbel_gompertz_lower_ci_theta)
gumbel_gompertz_lower_ci_rho <- rho(gumbel_gompertz_lower_ci_theta_copula)

gumbel_gompertz_g1 <- gumbel_gompertz_optim$par[1]
gumbel_gompertz_l1 <- gumbel_gompertz_optim$par[2]
gumbel_gompertz_g2 <- gumbel_gompertz_optim$par[3]
gumbel_gompertz_l2 <- gumbel_gompertz_optim$par[4]
gumbel_gompertz_theta <- gumbel_gompertz_optim$par[5]
gumbel_gompertz_copula <- gumbelCopula(gumbel_gompertz_theta)
gumbel_gompertz_rho <- rho(gumbel_gompertz_copula)

#AIC
gumbel_gompertz_estimated_parameters <- c(gumbel_gompertz_g1, gumbel_gompertz_l1, gumbel_gompertz_g2, gumbel_gompertz_l2, gumbel_gompertz_theta)
gumbel_gompertz_loglik_estimated_parameters <- gumbel_gompertz_loglik(gumbel_gompertz_estimated_parameters,X,Y,d1,d2)
k <- length(gumbel_gompertz_estimated_parameters)
gumbel_gompertz_aic <- -2*gumbel_gompertz_loglik_estimated_parameters+2*k

table_S4$gamma1[4] <- gumbel_gompertz_g1
table_S4$gamma1_lwci[4] <- gumbel_gompertz_lower_ci_g1
table_S4$gamma1_upci[4] <- gumbel_gompertz_upper_ci_g1
table_S4$lambda1[4] <- gumbel_gompertz_l1
table_S4$lambda1_lwci[4] <- gumbel_gompertz_lower_ci_l1
table_S4$lambda1_upci[4] <- gumbel_gompertz_upper_ci_l1
table_S4$gamma2[4] <- gumbel_gompertz_g2
table_S4$gamma2_lwci[4] <- gumbel_gompertz_lower_ci_g2
table_S4$gamma2_upci[4] <- gumbel_gompertz_upper_ci_g2
table_S4$lambda2[4] <- gumbel_gompertz_l2
table_S4$lambda2_lwci[4] <- gumbel_gompertz_lower_ci_l2
table_S4$lambda2_upci[4] <- gumbel_gompertz_upper_ci_l2
table_S4$rho[4] <- gumbel_gompertz_rho
table_S4$rho_lwci[4] <- gumbel_gompertz_lower_ci_rho
table_S4$rho_upci[4] <- gumbel_gompertz_upper_ci_rho
table_S4$AIC[4] <- gumbel_gompertz_aic


##Normal & Gompertz copula analysis##
normal_gompertz_optim <- optim(c(0.02,0.1,0.02,0.1,0.4), normal_gompertz_loglik, method="L-BFGS-B",
                               lower=c(0.01,0.01,0.01,0.01,0.1),upper=c(0.1,0.1,0.6,0.1,0.9), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)

normal_gompertz_fisher_info <- solve(-normal_gompertz_optim$hessian) 
normal_gompertz_se <- sqrt(diag(normal_gompertz_fisher_info)) 

#95% CI for hazard rates and association parameter
normal_gompertz_upper_ci_g1 <- normal_gompertz_optim$par[1]+1.96*normal_gompertz_se[1]
normal_gompertz_lower_ci_g1 <- normal_gompertz_optim$par[1]-1.96*normal_gompertz_se[1]
normal_gompertz_upper_ci_l1 <- normal_gompertz_optim$par[2]+1.96*normal_gompertz_se[2]
normal_gompertz_lower_ci_l1 <- normal_gompertz_optim$par[2]-1.96*normal_gompertz_se[2]
normal_gompertz_upper_ci_g2 <- normal_gompertz_optim$par[3]+1.96*normal_gompertz_se[3]
normal_gompertz_lower_ci_g2 <- normal_gompertz_optim$par[3]-1.96*normal_gompertz_se[3]
normal_gompertz_upper_ci_l2 <- normal_gompertz_optim$par[4]+1.96*normal_gompertz_se[4]
normal_gompertz_lower_ci_l2 <- normal_gompertz_optim$par[4]-1.96*normal_gompertz_se[4]
normal_gompertz_upper_ci_theta <- normal_gompertz_optim$par[5]+1.96*normal_gompertz_se[5]
normal_gompertz_lower_ci_theta <- normal_gompertz_optim$par[5]-1.96*normal_gompertz_se[5]

#95% CI for association parameter to rho
normal_gompertz_upper_ci_theta_copula <- normalCopula(normal_gompertz_upper_ci_theta) #convert CI to rho
normal_gompertz_upper_ci_rho <- rho(normal_gompertz_upper_ci_theta_copula)
normal_gompertz_lower_ci_theta_copula <- normalCopula(normal_gompertz_lower_ci_theta)
normal_gompertz_lower_ci_rho <- rho(normal_gompertz_lower_ci_theta_copula)

normal_gompertz_g1 <- normal_gompertz_optim$par[1]
normal_gompertz_l1 <- normal_gompertz_optim$par[2]
normal_gompertz_g2 <- normal_gompertz_optim$par[3]
normal_gompertz_l2 <- normal_gompertz_optim$par[4]
normal_gompertz_theta <- normal_gompertz_optim$par[5]
normal_gompertz_copula <- normalCopula(normal_gompertz_theta)
normal_gompertz_rho <- rho(normal_gompertz_copula)

#AIC
normal_gompertz_estimated_parameters <- c(normal_gompertz_g1, normal_gompertz_l1, normal_gompertz_g2, normal_gompertz_l2, normal_gompertz_theta)
normal_gompertz_loglik_estimated_parameters <- normal_gompertz_loglik(normal_gompertz_estimated_parameters,X,Y,d1,d2)
k <- length(normal_gompertz_estimated_parameters)
normal_gompertz_aic <- -2*normal_gompertz_loglik_estimated_parameters+2*k

table_S4$gamma1[1] <- normal_gompertz_g1
table_S4$gamma1_lwci[1] <- normal_gompertz_lower_ci_g1
table_S4$gamma1_upci[1] <- normal_gompertz_upper_ci_g1
table_S4$lambda1[1] <- normal_gompertz_l1
table_S4$lambda1_lwci[1] <- normal_gompertz_lower_ci_l1
table_S4$lambda1_upci[1] <- normal_gompertz_upper_ci_l1
table_S4$gamma2[1] <- normal_gompertz_g2
table_S4$gamma2_lwci[1] <- normal_gompertz_lower_ci_g2
table_S4$gamma2_upci[1] <- normal_gompertz_upper_ci_g2
table_S4$lambda2[1] <- normal_gompertz_l2
table_S4$lambda2_lwci[1] <- normal_gompertz_lower_ci_l2
table_S4$lambda2_upci[1] <- normal_gompertz_upper_ci_l2
table_S4$rho[1] <- normal_gompertz_rho
table_S4$rho_lwci[1] <- normal_gompertz_lower_ci_rho
table_S4$rho_upci[1] <- normal_gompertz_upper_ci_rho
table_S4$AIC[1] <- normal_gompertz_aic


##Contour plots##

########################
## Recreating Figure 2## 
########################

#Set up environment
par(mfrow=c(2,2))
y <- c(0,0.2,0.4,0.6,0.8,1)

#Normal copula contour plot
t1 <- iRho(normalCopula(),0.634) #theta
myMvd1 <- mvdc(normalCopula(t1),
               margins = c("unif", "unif"), paramMargins = list(list(min=0,max=1),list(min=0,max=1))) #Normal copula with uniform marginal distributions
contour(myMvd1, dMvdc,xlim=c(0,1), ylim=c(0,1),yaxt="n",xlab="Survival probability of SI switch", ylab="Survival probability of death from AIDS",nlevels=14,cex.axis=1.3, cex.lab=1.3)
title("Normal copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Clayton copula contour plot
t2 <- iRho(claytonCopula(),0.754) #theta
myMvd2 <- mvdc(copula = archmCopula(family = "clayton", param = t2),
               margins = c("unif", "unif"), paramMargins=list(list(min=0,max=1),list(min=0,max=1))) #Clayton copula with uniform marginal distributions
contour(myMvd2, dMvdc,xlim=c(0,1), ylim=c(0,1),yaxt="n",xlab="Survival probability of SI switch", ylab="Survival probability of death from AIDS",nlevels=20, cex.axis=1.3, cex.lab=1.3)
title("Clayton copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Frank copula contour plot
t3 <- iRho(frankCopula(),0.697) #theta
myMvd3 <- mvdc(copula = archmCopula(family = "frank", param = t3),
               margins = c("unif", "unif"), paramMargins=list(list(min=0,max=1),list(min=0,max=1))) #Frank copula with uniform marginal distributions
contour(myMvd3, dMvdc,xlim=c(0,1), ylim=c(0,1),yaxt="n",xlab="Survival probability of SI switch", ylab="Survival probability of death from AIDS",nlevels=20,cex.axis=1.3, cex.lab=1.3)
title("Frank copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Gumbel copula contour plot
t4 <- iRho(gumbelCopula(),0.504) #theta
myMvd4 <- mvdc(copula = archmCopula(family = "gumbel", param = t4),
               margins = c("unif", "unif"), paramMargins=list(list(min=0,max=1),list(min=0,max=1))) #Gumbel copula with uniform marginal distributions
contour(myMvd4, dMvdc,xlim=c(0,1), ylim=c(0,1),xlab="Survival probability of SI switch", ylab="Survival probability of death from AIDS",nlevels=20,cex.axis=1.3, cex.lab=1.3)
title("Gumbel copula")
axis(2, at=y,labels=y, cex.axis=1.3)



##########################
## Recreating Figure S2 ## 
##########################

#Set up environment
par(mfrow=c(2,2))
y=c(0,2,4,6,8,10,12)

#Normal copula contour
t <- iRho(normalCopula(),0.6338) #theta
l1 <- 0.072 #hazard rate for non-terminal event
l2 <- 0.071 #hazard rate for terminal event
myMvd2 <- mvdc(normalCopula(t),
               margins = c("exp", "exp"), paramMargins = list(lambda=l1, lambda=l2)) #copula with Exponential marginal distributions
contour(myMvd2, dMvdc,xlim=c(0,12.5), ylim=c(0,12.5), yaxt="n",xlab="Time to SI switch", ylab="Time to death from AIDS",nlevels=14, cex.axis=1.2, cex.lab=1.3)
title("Normal copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Clayton copula contour
t1<-iRho(claytonCopula(),0.754) #theta
l1<-0.073 #hazard rate for non-terminal event
l2<-0.073 #hazard rate for terminal event
myMvd3 <- mvdc(copula = archmCopula(family = "clayton", param = t1),
               margins = c("exp", "exp"), paramMargins = list(lambda=l1, lambda=l2))
contour(myMvd3, dMvdc,xlim=c(0,12.5), ylim=c(0,12.5), yaxt="n",xlab="Time to SI switch", ylab="Time to death from AIDS",nlevels=40, cex.axis=1.2, cex.lab=1.3)
title("Clayton copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Frank copula contour
t2<-iRho(frankCopula(),0.697) #theta
l1<-0.076 #hazard rate for non-terminal event
l2<-0.072 #hazard rate for terminal event
myMvd4 <- mvdc(copula = archmCopula(family = "frank", param = t2),
               margins = c("exp", "exp"), paramMargins = list(lambda=l1, lambda=l2))
contour(myMvd4, dMvdc,xlim=c(0,12.5), ylim=c(0,12.5), yaxt="n",xlab="Time to SI switch", ylab="Time to death from AIDS",nlevels=14, cex.axis=1.2, cex.lab=1.3)
title("Frank copula")
axis(2, at=y,labels=y, cex.axis=1.3)

#Gumbel copula contour
t3<-iRho(gumbelCopula(),0.504) #theta
l1<-0.072 #hazard rate for non-terminal event
l2<-0.069 #hazard rate for terminal event
myMvd5 <- mvdc(copula = archmCopula(family = "gumbel", param = t3),
               margins = c("exp", "exp"), paramMargins = list(lambda=l1, lambda=l2))
contour(myMvd5, dMvdc,xlim=c(0,12.5), yaxt="n",ylim=c(0,12.5), xlab="Time to SI switch", ylab="Time to death from AIDS", nlevels=14, cex.axis=1.2, cex.lab=1.3)
title("Gumbel copula")
axis(2, at=y,labels=y, cex.axis=1.3)



