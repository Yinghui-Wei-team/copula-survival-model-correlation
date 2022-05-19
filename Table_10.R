#######################################################################################
#
#   Filename    :	      Table 10.R    												  
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

#prepare table
exponential_MAE <- rep(NA, 12)
exponential_CP <- rep(NA, 12)
weibull_MAE <- rep(NA, 12)
weibull_CP <- rep(NA, 12)
gompertz_MAE <- rep(NA, 12)
gompertz_CP <- rep(NA, 12)
copula <- c(rep("Normal",3), rep("Clayton",3), rep("Frank",3), rep("Gumbel",3))
survival <- rep(c("Exponential", "Weibull", "Gompertz"),4)
table_10 <- data.frame(copula, survival,
                       exponential_MAE, exponential_CP,  
                       weibull_MAE, weibull_CP,  
                       gompertz_MAE, gompertz_CP)

n <- 329                           #number of individuals 
runs <- 1000                       #number of data sets




### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(4210097)
true_theta <- 2.55
true_rho <- rho(claytonCopula(true_theta)) 
true_lambda1 <- 0.073                   #hazard rate for SI switch
true_lambda2 <- 0.073                   #hazard rate for death

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.03,0.03,2), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),
                    upper=c(0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci to rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[3]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[4] <- bias
table_10$exponential_CP[4] <- coverage




### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(65322677)
true_theta <- 2.55
true_rho <- rho(claytonCopula(true_theta)) 
true_lambda1 <- 0.073                   #hazard rate of SI switch
true_lambda2 <- 0.073                   #hazard rate of death

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[4] <- bias
table_10$weibull_CP[4] <- coverage



### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(8666857)
true_theta <- 2.55
true_rho <- rho(claytonCopula(true_theta)) 
true_lambda1 <- 0.073                   #hazard rate of SI switch
true_lambda2 <- 0.073                   #hazard rate of death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.1,0.5,0.1,0.5,2), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.001,-0.1,0.001,0.01),
                    upper=c(0.1,0.5,0.1,0.5,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[4] <- bias
table_10$gompertz_CP[4] <- coverage



### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(347444377)
true_theta <- 1.72
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 1.13                   #Weibull parameters for SI switch
true_b1 <- 0.05
true_a2 <- 1.95                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.1,0.1,1), clayton_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[3]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[5] <- bias
table_10$exponential_CP[5] <- coverage


### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(64563475)
true_theta <- 1.72
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 1.13        #Weibull parameters for SI switch           
true_b1 <- 0.05
true_a2 <- 1.95        #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #bounds for optimisation
  a1_lw <- 0.01
  a1_up <- 1.7
  b1_lw <- 0.001
  b1_up <- 0.2
  a2_lw <- 0.01
  a2_up <- 2.3
  b2_lw <- 0.001
  b2_up <- 0.2
  t_lw <- 0.01
  t_up <- 10
  
  clayton_optim <- optim(c(true_a1, true_b1, true_a2, true_b2, true_theta), clayton_weibull_loglik_sim, method="L-BFGS-B",lower=c(a1_lw, b1_lw,a2_lw,b2_lw,t_lw),
                    upper=c(a1_up, b1_up,a2_up,b2_up,t_up), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[5] <- bias
table_10$weibull_CP[5] <- coverage



### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(65321127)
true_theta <- 1.72
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 1.13                   #Weibull parameters for SI Switch
true_b1 <- 0.05
true_a2 <- 1.95                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.07, 0.05, 0.21, 0.02, 1.57), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[5] <- bias
table_10$gompertz_CP[5] <- coverage



### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(39985577)
true_theta <- 1.57
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.07                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.21                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.1,0.1,1), clayton_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[3]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[6] <- bias
table_10$exponential_CP[6] <- coverage




### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(64999887)
true_theta <- 1.57
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.07                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.21                  #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(1.13,0.05,1.95,0.01,2), clayton_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),
                    upper=c(1.5,0.2,3,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  if (lower_ci_theta< -1){lower_ci_theta=-1}
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[6] <- bias
table_10$weibull_CP[6] <- coverage




### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(65111177)
true_theta <- 1.57
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.07                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.21                  #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- claytonCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  clayton_optim <- optim(c(0.07, 0.05, 0.2, 0.02,1.5), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.02,0.01,-0.01,0.001,0.01),
                    upper=c(2,0.2,3,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  c.cop.up <- claytonCopula(upper_ci_theta) 
  upper_ci_rho <- rho(c.cop.up)
  c.cop.low <- claytonCopula(lower_ci_theta)
  lower_ci_rho <- rho(c.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- clayton_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- clayton_optim$par[5]
  c.cop <- claytonCopula(est_t_loop)
  estimated_rho[i] <- rho(c.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[6] <- bias
table_10$gompertz_CP[6] <- coverage




### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(34799097)
true_theta <- 5.8
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.076                   #hazard rate for SI switch
true_lambda2 <- 0.072                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,1), frank_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,23.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- frank_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[3]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[7] <- bias
table_10$exponential_CP[7] <- coverage



### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(35990957)
true_theta <- 5.8
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.076                   #hazard rate for SI switch
true_lambda2 <- 0.072                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,2.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
 
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[7] <- bias
table_10$weibull_CP[7] <- coverage



### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(3666545)
true_theta <- 5.8
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.076                   #hazard rate for SI switch
true_lambda2 <- 0.072                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,0.01),
                    upper=c(0.1,0.2,0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[7] <- bias
table_10$gompertz_CP[7] <- coverage



### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(8066477)
true_theta <- 5.42
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 1.06                   #Weibull parameters for SI switch
true_b1 <- 0.06                   
true_a2 <- 1.90                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,1), frank_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,23.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- frank_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[3]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[8] <- bias
table_10$exponential_CP[8] <- coverage



### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(8009895)
true_theta <- 5.42
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 1.06                   #Weibull parameters for SI switch
true_b1 <- 0.06                   
true_a2 <- 1.90                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),
                    upper=c(2,0.2,3,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)

  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[8] <- bias
table_10$weibull_CP[8] <- coverage



### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(93838945)
true_theta <- 5.42
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 1.06                   #Weibull parameters for SI switch
true_b1 <- 0.06                   
true_a2 <- 1.90                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.001,-0.1,0.001,0.01),
                    upper=c(0.1,0.4,0.5,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[8] <- bias
table_10$gompertz_CP[8] <- coverage



### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(67221904)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.05                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05                   
true_gamma2 <- 0.22                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,1), frank_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,23.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error 
  
  upper_ci_theta <- frank_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- frank_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[3]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[9] <- bias
table_10$exponential_CP[9] <- coverage




### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(8800986)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.05                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05                   
true_gamma2 <- 0.22                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(1,0.1,2,0.1,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.1,0.001,0.1,0.001,0.01),
                    upper=c(2,0.2,3,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[9] <- bias
table_10$weibull_CP[9] <- coverage



### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(90093121)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.05                   #Gompertz parameters for SI switch
true_lambda1 <- 0.05                   
true_gamma2 <- 0.22                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- frankCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  frank_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.001,0.01),
                    upper=c(0.2,0.1,0.5,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-frank_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- frank_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- frank_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  f.cop.up <- frankCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  f.cop.low <- frankCopula(lower_ci_theta)
  lower_ci_rho <- rho(f.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- frank_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- frank_optim$par[5]
  f.cop <- frankCopula(est_t_loop)
  estimated_rho[i] <- rho(f.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[9] <- bias
table_10$gompertz_CP[9] <- coverage





### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(24552112)
true_theta <- 1.55
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.069                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if (lower_ci_theta <1) {lower_ci_theta=1}
  
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[3]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[10] <- bias
table_10$exponential_CP[10] <- coverage




### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(80098509)
true_theta <- 1.55
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.069                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim <- optim(c(0.2,0.2,0.2,0.2,2), gumbel_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,1.01),
                    upper=c(2,1,3,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[10] <- bias
table_10$weibull_CP[10] <- coverage




### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(982586609)
true_theta <- 1.55
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.069                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim <- optim(c(0.01,0.04,0.01,0.03,1.17), gumbel_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,1.01),
                    upper=c(0.2,0.2,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[10] <- bias
table_10$gompertz_CP[10] <- coverage



### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(287779545)
true_theta <- 1.40
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.91          #Weibull parameters for SI switch         
true_b1 <- 0.08
true_a2 <- 1.93          #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of thetagumbel_optim
  estimated_theta[i] <- gumbel_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[3]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[11] <- bias
table_10$exponential_CP[11] <- coverage



### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(1997536)
true_theta <- 1.40
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.91                   #Weibull parameters for SI switch
true_b1 <- 0.08
true_a2 <- 1.93                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim  <- optim(c(0.91,0.08,1.93,0.01,2), gumbel_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,1.01),
                    upper=c(2,1,3,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[11] <- bias
table_10$weibull_CP[11] <- coverage




### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(79732115)
true_theta <- 1.40
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.91              #Weibull parameters for SI switch     
true_b1 <- 0.08
true_a2 <- 1.93              #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim  <- optim(c(0.07,0.05,0.21,0.02,1.17), gumbel_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,1.01),
                    upper=c(0.2,0.2,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error 
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[11] <- bias
table_10$gompertz_CP[11] <- coverage



### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(7009808)
true_theta <- 1.44
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.03                   #Gompertz parameters for SI switch
true_lambda1 <- 0.06
true_gamma2 <- 0.23                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[3]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[12] <- bias
table_10$exponential_CP[12] <- coverage




### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(2406645)
true_theta <- 1.44
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.03                   #Gompertz parameters for SI switch
true_lambda1 <- 0.06
true_gamma2 <- 0.23                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  a1_lw <- 0.2
  a1_up <- 1.9
  b1_lw <- 0.01
  b1_up <- 0.3
  a2_lw <- 0.2
  a2_up <- 3
  b2_lw <- 0.00001
  b2_up <- 0.7
  t_lw <- 1.01
  t_up <- 10
  
  gumbel_optim <- optim(c(0.6,0.06,0.84,0.04,1.26), gumbel_weibull_loglik_sim, method="L-BFGS-B",lower=c(a1_lw, b1_lw, a2_lw, b2_lw, t_lw),
                    upper=c(a1_up, b1_up, a2_up, b2_up, t_up), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[12] <- bias
table_10$weibull_CP[12] <- coverage



### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(50490495)
true_theta <- 1.44
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.03                   #Gompertz parameters for SI switch
true_lambda1 <- 0.06
true_gamma2 <- 0.23                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- gumbelCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  gumbel_optim  <- optim(c(0.01,0.04,0.01,0.03,1.17), gumbel_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.001,1.01),
                    upper=c(0.2,0.2,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  
  #convert theta ci -> rho ci
  g.cop.up <- gumbelCopula(upper_ci_theta) 
  upper_ci_rho <- rho(g.cop.up)
  g.cop.low <- gumbelCopula(lower_ci_theta)
  lower_ci_rho <- rho(g.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- gumbel_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- gumbel_optim$par[5]
  g.cop <- gumbelCopula(est_t_loop)
  estimated_rho[i] <- rho(g.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100

table_10$gompertz_MAE[12] <- bias
table_10$gompertz_CP[12] <- coverage





### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(2900868)
true_theta <- 0.65
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.071                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.01,0.01,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.1,0.1,0.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- normal_optim$par[3]-1.96*se[3] #lower ci for theta
  
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[3]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[1] <- bias
table_10$exponential_CP[1] <- coverage



### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(80006065)
true_theta <- 0.65
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.071                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.1),upper=c(2,0.2,3,0.2,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[1] <- bias
table_10$weibull_CP[1] <- coverage



### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(90080098)
true_theta <- 0.65
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.072                   #hazard rate for SI switch
true_lambda2 <- 0.071                   #hazard rate for death from AIDS

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- -log(U1)/true_lambda1 #time to non-terminal event, using inverse Exponential survival function 
  T2 <- -log(V1)/true_lambda2 #time to terminal event, using inverse Exponential survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.01,0.01,-0.03,0.01,0.01),upper=c(0.2,0.1,0.9,0.1,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
 
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[1] <- bias
table_10$gompertz_CP[1] <- coverage




### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(89987642)
true_theta <- 0.53
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 1.21                    #Weibull parameters for SI switch
true_b1 <- 0.04
true_a2 <- 1.96                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.2,0.2,0.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error 
  
  upper_ci_theta <- normal_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- normal_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(f.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[3]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[2] <- bias
table_10$exponential_CP[2] <- coverage



### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(43326765)
true_theta <- 0.53
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 1.21                    #Weibull parameters for SI switch
true_b1 <- 0.04
true_a2 <- 1.96                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim  <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.1),upper=c(2,0.1,3,0.1,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
 
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[2] <- bias
table_10$weibull_CP[2] <- coverage




### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(99868799)
true_theta <- 0.53
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 1.21                   #Weibull parameters for SI switch
true_b1 <- 0.04
true_a2 <- 1.96                   #Weibull parameters for death from AIDS
true_b2 <- 0.01

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- (-log(U1)/true_b1)^(1/true_a1) #time to non-terminal event, using inverse Weibull survival function 
  T2 <- (-log(V1)/true_b2)^(1/true_a2) #time to terminal event, using inverse Weibull survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim  <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.01,0.01,-0.03,0.01,0.01),upper=c(0.2,0.1,0.9,0.1,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[2] <- bias
table_10$gompertz_CP[2] <- coverage




### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(61121885)
true_theta <- 0.57
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.05                    #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.23                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.1,0.1,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.1,0.1,0.9),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- normal_optim$par[3]-1.96*se[3] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[3]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[3]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$exponential_MAE[3] <- bias
table_10$exponential_CP[3] <- coverage




### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(588002709)
true_theta <- 0.59
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.05                    #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.23                   #Gompertz parameter for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.1),upper=c(2,0.1,3,0.1,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lower_ci_rho <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$weibull_MAE[3] <- bias
table_10$weibull_CP[3] <- coverage




### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(8009668)
true_theta <- 0.59
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.05                    #Gompertz parameters for SI switch
true_lambda1 <- 0.05
true_gamma2 <- 0.23                   #Gompertz parameters for death from AIDS
true_lambda2 <- 0.02

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <- rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:(n)){   #loop to simulate U an V
    m=1                  
    u1 <- runif(m,0,1) #uniformly transformed time to non-terminal event
    fc <- normalCopula(true_theta, dim=2) #conditional distribution method
    uv <- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) 
    u <- uv[,1]  #uniformly transformed time to non-terminal event
    v <- uv[,2]  #uniformly transformed time to terminal event
    U1[k]=u      #add to u and v vectors on the outside of loop
    V1[k]=v
  }
  T1 <- 1/true_gamma1 *log (1-true_gamma1/true_lambda1 *log(U1)) #time to non-terminal event, using inverse Gompertz survival function 
  T2 <- 1/true_gamma2 *log (1-true_gamma2/true_lambda2 *log(V1)) #time to terminal event, using inverse Gompertz survival function 
  C <- runif(n,0,50) #censoring time
  X <- pmin(T1,T2,C) #observed time to non-terminal event
  Y <- pmin(T2, C)  #observed time to terminal event
  d1 <- ifelse(T1<=Y,1,0) #event indicator for non-terminal event
  d2 <- ifelse(T2<=C,1,0) #event indicator for terminal event
  df <- data.frame(X, Y, d1, d2) 
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  #optimise
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.01,0.01,-0.03,0.01,0.01),upper=c(0.2,0.1,0.9,0.1,0.9), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-normal_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- normal_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- normal_optim$par[5]-1.96*se[5] #lower ci for theta
  
  #convert theta ci -> rho ci
  n.cop.up <- normalCopula(upper_ci_theta) 
  upper_ci_rho <- rho(n.cop.up)
  n.cop.low <- normalCopula(lower_ci_theta)
  lclower_ci_rhoi3 <- rho(n.cop.low)
  
  #coverage for rho
  if(true_rho <= upper_ci_rho && true_rho >= lower_ci_rho) {counter_rho = counter_rho+1}
  lower_ci[i] <- lower_ci_rho
  upper_ci[i] <- upper_ci_rho
  
  #save estimates of theta 
  estimated_theta[i] <- normal_optim$par[5]
  
  #save estimates of rho
  est_t_loop <- normal_optim$par[5]
  n.cop <- normalCopula(est_t_loop)
  estimated_rho[i] <- rho(n.cop)
  
  print(i)
}

#Results
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_10$gompertz_MAE[3] <- bias
table_10$gompertz_CP[3] <- coverage
