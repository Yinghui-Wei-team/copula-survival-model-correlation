#######################################################################################
#
#   Filename    :	      Table 7.R    												  
#                                                                                     
#   Project     :       Article "Estimating the correlation between semi-competing risk survival endpoints"                                                             
#   Authors     :       L Sorrell, Y Wei, M Wojtys and P Rowe                                                               
#   Date        :       01/06/2021
#																				  
#   R Version   :       R-3.6.1                                                              
#
#   Required R packages :  copula, mvtnorm
#
########################################################################################
library(copula) 
library(mvtnorm)
source("Functions.R")   #Likelihood functions

#prepare table:
lambda1_MAE <- rep(NA, 17)
lambda1_CP <- rep(NA, 17)
lambda2_MAE <- rep(NA, 17)
lambda2_CP <- rep(NA, 17)
rho_MAE <- rep(NA, 17)
rho_CP <- rep(NA, 17)
Normal_perc <- rep(NA, 17)
Clayton_perc <- rep(NA, 17)
Frank_perc <- rep(NA, 17)
Gumbel_perc <- rep(NA, 17)
rho <- c(0, rep(0.25,4), rep(0.5,4), rep(0.75, 4), rep(0.9,4))
copula <- c("Independence", rep(c("Normal", "Clayton", "Frank", "Gumbel"),4))
table_7 <- data.frame(rho, copula, 
                      lambda1_MAE, lambda1_CP, 
                      lambda2_MAE, lambda2_CP, 
                      rho_MAE, rho_CP, 
                      Normal_perc, Clayton_perc, Frank_perc, Gumbel_perc)

true_lambda1 <- 0.05514            #hazard rate for SI switch
true_lambda2 <- 0.0725             #hazard rate for death from AIDS
n <- 329                           #number of individuals 
runs <- 1000                       #number of data sets



##Independence copula with rho=0##
set.seed(2947733)
true_theta <- 0
copula_normal <- normalCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_normal)

#Save & reset for reporting 
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  normal_copula <- normalCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = normal_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[1] <- bias_lambda1
table_7$lambda1_CP[1] <- coverage_lambda1
table_7$lambda2_MAE[1] <- bias_lambda2
table_7$lambda2_CP[1] <- coverage_lambda2
table_7$rho_MAE[1] <- bias_rho
table_7$rho_CP[1] <- coverage_rho
table_7$Normal_perc[1] <- percent_normal
table_7$Clayton_perc[1] <- percent_clayton
table_7$Frank_perc[1] <- percent_frank
table_7$Gumbel_perc[1] <- percent_gumbel






##Normal copula with rho=0.25##
set.seed(92736611)
true_theta <- 0.26
copula_normal <- normalCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_normal)

#Save & reset for reporting 
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  normal_copula <- normalCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = normal_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[2] <- bias_lambda1
table_7$lambda1_CP[2] <- coverage_lambda1
table_7$lambda2_MAE[2] <- bias_lambda2
table_7$lambda2_CP[2] <- coverage_lambda2
table_7$rho_MAE[2] <- bias_rho
table_7$rho_CP[2] <- coverage_rho
table_7$Normal_perc[2] <- percent_normal
table_7$Clayton_perc[2] <- percent_clayton
table_7$Frank_perc[2] <- percent_frank
table_7$Gumbel_perc[2] <- percent_gumbel





##Clayton copula with rho=0.25##
set.seed(93726289)
true_theta <- 0.4
copula_clayton <- claytonCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_clayton)

#Save & reset for reporting 
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  clayton_copula <- claytonCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = clayton_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[3] <- bias_lambda1
table_7$lambda1_CP[3] <- coverage_lambda1
table_7$lambda2_MAE[3] <- bias_lambda2
table_7$lambda2_CP[3] <- coverage_lambda2
table_7$rho_MAE[3] <- bias_rho
table_7$rho_CP[3] <- coverage_rho
table_7$Normal_perc[3] <- percent_normal
table_7$Clayton_perc[3] <- percent_clayton
table_7$Frank_perc[3] <- percent_frank
table_7$Gumbel_perc[3] <- percent_gumbel







##Frank copula with rho=0.25##
set.seed(928473381)
true_theta <- 1.54
copula_frank <- frankCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_frank)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  frank_copula <- frankCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = frank_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1

  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[4] <- bias_lambda1 
table_7$lambda1_CP[4] <- coverage_lambda1 
table_7$lambda2_MAE[4] <- bias_lambda2
table_7$lambda2_CP[4] <- coverage_lambda2
table_7$rho_MAE[4] <- bias_rho
table_7$rho_CP[4] <- coverage_rho
table_7$Normal_perc[4] <- percent_normal
table_7$Clayton_perc[4] <- percent_clayton
table_7$Frank_perc[4] <- percent_frank
table_7$Gumbel_perc[4] <- percent_gumbel





##Gumbel copula with rho=0.25##
set.seed(679137113)
true_theta <- 1.2
copula_gumbel <- gumbelCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_gumbel)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  gumbel_copula <- gumbelCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = gumbel_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1

  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[5] <- bias_lambda1
table_7$lambda1_CP[5] <- coverage_lambda1
table_7$lambda2_MAE[5] <- bias_lambda2
table_7$lambda2_CP[5] <- coverage_lambda2
table_7$rho_MAE[5] <- bias_rho
table_7$rho_CP[5] <- coverage_rho
table_7$Normal_perc[5] <- percent_normal
table_7$Clayton_perc[5] <- percent_clayton
table_7$Frank_perc[5] <- percent_frank
table_7$Gumbel_perc[5] <- percent_gumbel




##Normal copula with rho=0.5##
set.seed(912873912)
true_theta <- 0.52
copula_normal <- normalCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_normal)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  normal_copula <- normalCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = normal_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[6] <- bias_lambda1
table_7$lambda1_CP[6] <- coverage_lambda1
table_7$lambda2_MAE[6] <- bias_lambda2
table_7$lambda2_CP[6] <- coverage_lambda2
table_7$rho_MAE[6] <- bias_rho
table_7$rho_CP[6] <- coverage_rho
table_7$Normal_perc[6] <- percent_normal
table_7$Clayton_perc[6] <- percent_clayton
table_7$Frank_perc[6] <- percent_frank
table_7$Gumbel_perc[6] <- percent_gumbel





##Clayton copula with rho=0.5##
set.seed(12897686)
true_theta <- 3.45
copula_clayton <- claytonCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_clayton)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  clayton_copula <- claytonCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = clayton_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[7] <- bias_lambda1
table_7$lambda1_CP[7] <- coverage_lambda1
table_7$lambda2_MAE[7] <- bias_lambda2
table_7$lambda2_CP[7] <- coverage_lambda2
table_7$rho_MAE[7] <- bias_rho
table_7$rho_CP[7] <- coverage_rho
table_7$Normal_perc[7] <- percent_normal
table_7$Clayton_perc[7] <- percent_clayton
table_7$Frank_perc[7] <- percent_frank
table_7$Gumbel_perc[7] <- percent_gumbel







##Frank copula with rho=0.5##
set.seed(21937297)
true_theta <- 1.08
copula_frank <- frankCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_frank)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  frank_copula <- frankCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = frank_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[8] <- bias_lambda1
table_7$lambda1_CP[8] <- coverage_lambda1
table_7$lambda2_MAE[8] <- bias_lambda2
table_7$lambda2_CP[8] <- coverage_lambda2
table_7$rho_MAE[8] <- bias_rho
table_7$rho_CP[8] <- coverage_rho
table_7$Normal_perc[8] <- percent_normal
table_7$Clayton_perc[8] <- percent_clayton
table_7$Frank_perc[8] <- percent_frank
table_7$Gumbel_perc[8] <- percent_gumbel





##Gumbel copula with rho=0.5##
set.seed(129873017)
true_theta <- 1.54
copula_gumbel <- gumbelCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_gumbel)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  gumbel_copula <- gumbelCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = gumbel_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[9] <- bias_lambda1
table_7$lambda1_CP[9] <- coverage_lambda1
table_7$lambda2_MAE[9] <- bias_lambda2
table_7$lambda2_CP[9] <- coverage_lambda2
table_7$rho_MAE[9] <- bias_rho
table_7$rho_CP[9] <- coverage_rho
table_7$Normal_perc[9] <- percent_normal
table_7$Clayton_perc[9] <- percent_clayton
table_7$Frank_perc[9] <- percent_frank
table_7$Gumbel_perc[9] <- percent_gumbel






##Normal copula with rho=0.75##
set.seed(123481394)
true_theta <- 0.77
copula_normal <- normalCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_normal)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  normal_copula <- normalCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = normal_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim  <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99),X=df$X,Y=df$Y,d1=df$d1,d2=df$d2 ,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[10] <- bias_lambda1
table_7$lambda1_CP[10] <- coverage_lambda1
table_7$lambda2_MAE[10] <- bias_lambda2
table_7$lambda2_CP[10] <- coverage_lambda2
table_7$rho_MAE[10] <- bias_rho
table_7$rho_CP[10] <- coverage_rho
table_7$Normal_perc[10] <- percent_normal
table_7$Clayton_perc[10] <- percent_clayton
table_7$Frank_perc[10] <- percent_frank
table_7$Gumbel_perc[10] <- percent_gumbel





##Clayton copula with rho=0.75##
set.seed(91296396)
true_theta <- 2.58
copula_clayton <- claytonCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_clayton)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  clayton_copula <- claytonCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = clayton_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[11] <- bias_lambda1
table_7$lambda1_CP[11] <- coverage_lambda1
table_7$lambda2_MAE[11] <- bias_lambda2
table_7$lambda2_CP[11] <- coverage_lambda2
table_7$rho_MAE[11] <- bias_rho
table_7$rho_CP[11] <- coverage_rho
table_7$Normal_perc[11] <- percent_normal
table_7$Clayton_perc[11] <- percent_clayton
table_7$Frank_perc[11] <- percent_frank
table_7$Gumbel_perc[11] <- percent_gumbel






##Frank copula with rho=0.75##
set.seed(182739883)
true_theta <- 6.74
copula_frank <- frankCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_frank)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  frank_copula <- frankCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = frank_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[12] <- bias_lambda1
table_7$lambda1_CP[12] <- coverage_lambda1
table_7$lambda2_MAE[12] <- bias_lambda2
table_7$lambda2_CP[12] <- coverage_lambda2
table_7$rho_MAE[12] <- bias_rho
table_7$rho_CP[12] <- coverage_rho
table_7$Normal_perc[12] <- percent_normal
table_7$Clayton_perc[12] <- percent_clayton
table_7$Frank_perc[12] <- percent_frank
table_7$Gumbel_perc[12] <- percent_gumbel







##Gumbel copula with rho=0.75##
set.seed(238947938)
true_theta <- 2.29
copula_gumbel <- gumbelCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_gumbel)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  gumbel_copula <- gumbelCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = gumbel_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[13] <- bias_lambda1
table_7$lambda1_CP[13] <- coverage_lambda1
table_7$lambda2_MAE[13] <- bias_lambda2
table_7$lambda2_CP[13] <- coverage_lambda2
table_7$rho_MAE[13] <- bias_rho
table_7$rho_CP[13] <- coverage_rho
table_7$Normal_perc[13] <- percent_normal
table_7$Clayton_perc[13] <- percent_clayton
table_7$Frank_perc[13] <- percent_frank
table_7$Gumbel_perc[13] <- percent_gumbel






##Normal copula with rho=0.9##
set.seed(969976)
true_theta <- 0.91
copula_normal <- normalCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_normal)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  normal_copula <- normalCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = normal_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[14] <- bias_lambda1
table_7$lambda1_CP[14] <- coverage_lambda1
table_7$lambda2_MAE[14] <- bias_lambda2
table_7$lambda2_CP[14] <- coverage_lambda2
table_7$rho_MAE[14] <- bias_rho
table_7$rho_CP[14] <- coverage_rho
table_7$Normal_perc[14] <- percent_normal
table_7$Clayton_perc[14] <- percent_clayton
table_7$Frank_perc[14] <- percent_frank
table_7$Gumbel_perc[14] <- percent_gumbel







##Clayton copula with rho=0.9##
set.seed(645246)
true_theta <- 5.58
copula_clayton <- claytonCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_clayton)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  clayton_copula <- claytonCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = clayton_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[15] <- bias_lambda1
table_7$lambda1_CP[15] <- coverage_lambda1
table_7$lambda2_MAE[15] <- bias_lambda2
table_7$lambda2_CP[15] <- coverage_lambda2
table_7$rho_MAE[15] <- bias_rho
table_7$rho_CP[15] <- coverage_rho
table_7$Normal_perc[15] <- percent_normal
table_7$Clayton_perc[15] <- percent_clayton
table_7$Frank_perc[15] <- percent_frank
table_7$Gumbel_perc[15] <- percent_gumbel






##Frank copula with rho=0.9##
set.seed(72112)
true_theta <- 12.25
copula_frank <- frankCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_frank)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  frank_copula <- frankCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = frank_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)

#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[16] <- bias_lambda1
table_7$lambda1_CP[16] <- coverage_lambda1
table_7$lambda2_MAE[16] <- bias_lambda2
table_7$lambda2_CP[16] <- coverage_lambda2
table_7$rho_MAE[16] <- bias_rho
table_7$rho_CP[16] <- coverage_rho
table_7$Normal_perc[16] <- percent_normal
table_7$Clayton_perc[16] <- percent_clayton
table_7$Frank_perc[16] <- percent_frank
table_7$Gumbel_perc[16] <- percent_gumbel







##Gumbel copula with rho=0.9##
set.seed(85329)
true_theta <- 3.74
copula_gumbel <- gumbelCopula(true_theta) #convert theta to rho
true_rho <- rho(copula_gumbel)

#save for reporting
counter_lambda1 = 0
counter_lambda2 = 0
counter_rho = 0
counter_normal = 0
counter_clayton = 0 
counter_frank = 0
counter_gumbel = 0
counternan = 0
counterless = 0
estimated_lambda1 <- rep(0,runs)
estimated_lambda2 <- rep(0,runs)
estimated_rho <- rep(0,runs)
true_lambda1_vec <- rep(true_lambda1, runs)
true_lambda2_vec <- rep(true_lambda2, runs)
true_rho_vec <- rep(true_rho, runs)

for(i in 1:(runs)){
  
  u1 <- runif(n,0,1) #uniformly transformed time to non-terminal event
  gumbel_copula <- gumbelCopula(true_theta, dim=2) 
  uv <- cCopula(cbind(u1, runif(n)), copula = gumbel_copula, inverse = TRUE) #conditional distribution method
  u <- uv[,1]  #uniformly transformed time to non-terminal event
  v <- uv[,2]  #uniformly transformed time to terminal event
  T1 <- -log(u)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(v)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,30) #Censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #Optimise with all copula models
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,0.99), X=df$X,Y=df$Y,d1=df$d1,d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  normal_estimated_params <- c(normal_optim$par[1], normal_optim$par[2],normal_optim$par[3]) 
  
  clayton_optim <- optim(c(0.1,0.1,0.4), clayton_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  clayton_estimated_params <- c(clayton_optim$par[1], clayton_optim$par[2],clayton_optim$par[3]) 
  
  frank_optim <- optim(c(0.1,0.1,0.4), frank_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  frank_estimated_params <- c(frank_optim$par[1], frank_optim$par[2],frank_optim$par[3]) 
  
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.01,0.001,1.005),upper=c(0.2,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  gumbel_estimated_params <- c(gumbel_optim$par[1], gumbel_optim$par[2],gumbel_optim$par[3]) 
  
  #AIC for each copula model
  loglik_gumbel <- gumbel_loglik(gumbel_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Gumbel copula model
  loglik_frank <- frank_loglik(frank_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)   #max log-likelihood for Frank copula model
  loglik_clayton <- clayton_loglik(clayton_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)              #max log-likelihood for Clayton copula model
  loglik_normal <- normal_loglik(normal_estimated_params, X=df$X,Y=df$Y,d1=df$d1,d2=df$d2)  #max log-likelihood for Normal copula model 
  
  k <- 3  #3 parameters: lambda1,lambda2,theta  
  aic_gumbel <- -2*loglik_gumbel + 2*k #AIC for Gumbel copula model 
  aic_frank <- -2*loglik_frank + 2*k #AIC for Frank copula model
  aic_clayton <- -2*loglik_clayton + 2*k #AIC for Clayton copula model
  aic_normal <- -2*loglik_normal + 2*k #AIC for Normal copula model
  
  aics <- c(aic_normal, aic_clayton, aic_frank, aic_gumbel)  
  cops <- c("Normal", "Clayton", "Frank","Gumbel") 
  index <- which.min(aics) #Index copula model with minimum AIC     
  print(cops[index]) #Print which copula                                              
  
  #Continue with chosen copula
  if (index==1) #Continue with Normal copula model
  { 
    fisher_info <- solve(-normal_optim$hessian)  
    lambda1_est <- normal_optim$par[1]
    lambda2_est <- normal_optim$par[2]
    theta_est <- normal_optim$par[3] 
  } else if (index==2) #Continue with Clayton copula model
  { 
    fisher_info <- solve(-clayton_optim$hessian)  
    lambda1_est <- clayton_optim$par[1]
    lambda2_est <- clayton_optim$par[2]
    theta_est <- clayton_optim$par[3]
  } else if (index==3) #Continue with Frank copula model
  { 
    fisher_info <- solve(-frank_optim$hessian)  
    lambda1_est <- frank_optim$par[1]
    lambda2_est <- frank_optim$par[2]
    theta_est <- frank_optim$par[3]
  } else  { #Continue with Gumbel copula model
    fisher_info <- solve(-gumbel_optim$hessian)  
    lambda1_est <- gumbel_optim$par[1]
    lambda2_est <- gumbel_optim$par[2]
    theta_est <- gumbel_optim$par[3]
  }
  
  se <- sqrt(diag(fisher_info))                  #standard error
  
  #95% CI for lambda1, lambda2 and theta
  upperci1 <- lambda1_est + 1.96*se[1]        #upper ci for lambda1
  lowerci1 <- lambda1_est - 1.96*se[1]        #lower ci for lambda1
  upperci2 <- lambda2_est + 1.96*se[2]        #upper ci for lambda2
  lowerci2 <- lambda2_est - 1.96*se[2]        #lower ci for lambda2
  upperci3 <- theta_est + 1.96*se[3]          #upper ci for theta
  lowerci3 <- theta_est - 1.96*se[3]          #lower ci for theta
  
  #Diagonal of hessian negative - count nans
  if (lowerci3=="NaN") {counternan=counternan+1}  
  if (lowerci3=="NaN") {next}  
  
  #Boundries of Clayton (below minimum of copula package bounds)
  if (index==2){
    if (lowerci3 < -1) {counterless=counterless+1}
    if (lowerci3 < -1) {lowerci3=-1} 
  }
  #Boundries of Gumbel (below minimum of copula package bounds)
  if (index==4){
    if (lowerci3 < 1) {counterless=counterless+1}
    if (lowerci3 < 1) {lowerci3=1}
  }
  
  #Convert theta to rho using copula package for chosen copula
  if (index==1) 
  {
    t.up <- normalCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- normalCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- normalCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==2)
  {
    t.up <- claytonCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- claytonCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- claytonCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else if (index==3)
  {
    t.up <- frankCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- frankCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- frankCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  } else 
  {
    t.up <- gumbelCopula(upperci3)          #convert upper CI for theta to rho
    rho_up <- rho(t.up) 
    
    t.low <- gumbelCopula(lowerci3)          #convert upper CI for theta to rho
    rho_low <- rho(t.low) 
    
    t.est <- gumbelCopula(theta_est)          #convert upper CI for theta to rho
    rho_est <- rho(t.est) 
  }
  
  #counting chosen copula
  if (index==1) {counter_normal = counter_normal+1}
  if (index==2) {counter_clayton = counter_clayton+1}
  if (index==3) {counter_frank = counter_frank+1}
  if (index==4) {counter_gumbel = counter_gumbel+1}
  
  #save estimates for biases
  estimated_lambda1[i] <- lambda1_est
  estimated_lambda2[i] <- lambda2_est
  estimated_rho[i] <- rho_est
  
  #Add to counter if true value is within confidence interval
  if(true_lambda1 <= upperci1 && true_lambda1 >= lowerci1) {counter_lambda1 = counter_lambda1 + 1}
  if(true_lambda2 <= upperci2 && true_lambda2 >= lowerci2) {counter_lambda2 = counter_lambda2 + 1}
  if(true_rho <= rho_up && true_rho >= rho_low)      {counter_rho = counter_rho + 1}
  
  print(i)
}

#correct copula
percent_normal <- (counter_normal / runs) * 100
percent_clayton <- (counter_clayton / runs) * 100
percent_frank <- (counter_frank / runs) * 100
percent_gumbel <- (counter_gumbel / runs) * 100

#mean absolute bias
bias_lambda1 <- mean(abs(estimated_lambda1-true_lambda1_vec))
bias_lambda2 <- mean(abs(estimated_lambda2-true_lambda2_vec))
bias_rho <- mean(abs(estimated_rho-true_rho_vec))

#mse
mse_lambda1 <- mean((estimated_lambda1-true_lambda1_vec)^2)
mse_lambda2 <-  mean((estimated_lambda2-true_lambda2_vec)^2)
mse_rho <- mean((estimated_rho-true_rho_vec)^2)
#coverage
coverage_lambda1 <-  (counter_lambda1 / runs) * 100
coverage_lambda2 <- (counter_lambda2 / runs) * 100
coverage_rho <- (counter_rho / runs) * 100

#save
table_7$lambda1_MAE[17] <- bias_lambda1
table_7$lambda1_CP[17] <- coverage_lambda1
table_7$lambda2_MAE[17] <- bias_lambda2
table_7$lambda2_CP[17] <- coverage_lambda2
table_7$rho_MAE[17] <- bias_rho
table_7$rho_CP[17] <- coverage_rho
table_7$Normal_perc[17] <- percent_normal
table_7$Clayton_perc[17] <- percent_clayton
table_7$Frank_perc[17] <- percent_frank
table_7$Gumbel_perc[17] <- percent_gumbel