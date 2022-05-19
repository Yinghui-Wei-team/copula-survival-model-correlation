#######################################################################################
#
#   Filename    :	      Table 9.R    												  
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

#prepare table
exponential_MAE <- rep(NA, 12)
exponential_CP <- rep(NA, 12)
weibull_MAE <- rep(NA, 12)
weibull_CP <- rep(NA, 12)
gompertz_MAE <- rep(NA, 12)
gompertz_CP <- rep(NA, 12)
copula <- c(rep("Normal",3), rep("Clayton",3), rep("Frank",3), rep("Gumbel",3))
survival <- rep(c("Exponential", "Weibull", "Gompertz"),4)
table_9 <- data.frame(copula, survival,
                      exponential_MAE, exponential_CP,  
                      weibull_MAE, weibull_CP,
                      gompertz_MAE, gompertz_CP)

n <- 1199                        #number of individuals 
runs <- 1000                     #number of data sets


### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(34775677) 
true_theta <- 2.3 #True theta when generating data from Exponential distribution
true_rho <- rho(claytonCopula(true_theta))
true_lambda1 <- 0.029                   #kidney failure
true_lambda2 <- 0.030                   #death

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


table_9$exponential_MAE[4] <- bias
table_9$exponential_CP[4] <- coverage




### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(64445677)
true_theta <- 2.3
true_rho <- rho(claytonCopula(true_theta))
true_lambda1 <- 0.029                   #kidney failure
true_lambda2 <- 0.030                   #death

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
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[4] <- bias
table_9$weibull_CP[4] <- coverage



### Clayton copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(89129087)
true_theta <- 2.3
true_rho <- rho(claytonCopula(true_theta))
true_lambda1 <- 0.029                   #kidney failure
true_lambda2 <- 0.030                   #death

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
  clayton_optim <- optim(c(0.1,0.5,0.1,0.5,2), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.001,-0.1,0.001,0.01),
                    upper=c(0.1,0.5,0.1,0.5,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[4] <- bias
table_9$gompertz_CP[4] <- coverage





### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(347753377)
true_theta <- 1.96
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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


table_9$exponential_MAE[5] <- bias
table_9$exponential_CP[5] <- coverage




### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(64563477)
true_theta <- 1.96
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[5] <- bias
table_9$weibull_CP[5] <- coverage





### Clayton copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(65530977)
true_theta <- 1.96
true_rho <- rho(claytonCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[5] <- bias
table_9$gompertz_CP[5] <- coverage





### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(345543377)
true_theta <- 2.15
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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


table_9$exponential_MAE[6] <- bias
table_9$exponential_CP[6] <- coverage





### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(6400477)
true_theta <- 1.96
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[6] <- bias
table_9$weibull_CP[6] <- coverage





### Clayton copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(34554336)
true_theta <- 2.15
true_rho <- rho(claytonCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  clayton_optim <- optim(c(0.2,0.2,0.2,0.2,2), clayton_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,0.01),
                    upper=c(1.5,0.2,1.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[6] <- bias
table_9$gompertz_CP[6] <- coverage





### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(347677)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.028                   #kidney failure
true_lambda2 <- 0.030                   #death

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


table_9$exponential_MAE <- bias
table_9$exponential_CP[7] <- coverage




### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(35564557)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.028                   #kidney failure
true_lambda2 <- 0.030                   #death

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.1,0.2,1.1,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[7] <- bias
table_9$weibull_CP[7] <- coverage



### Frank copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(3666545)
true_theta <- 4.27
true_rho <- rho(frankCopula(true_theta))
true_lambda1 <- 0.028                   #kidney failure
true_lambda2 <- 0.030                   #death

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
                    upper=c(0.1,0.1,0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[7] <- bias
table_9$gompertz_CP[7] <- coverage





### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(347677)
true_theta <- 3.82
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06                   
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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


table_9$exponential_MAE[8] <- bias
table_9$exponential_CP[8] <- coverage




### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(35549995)
true_theta <- 3.82
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06                   
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1,0.2,1,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[8] <- bias
table_9$weibull_CP[8] <- coverage





### Frank copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(45577745)
true_theta <- 3.82
true_rho <- rho(frankCopula(true_theta))
true_a1 <- 0.63                   #kidney failure
true_b1 <- 0.06                   
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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
  frank_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,0.01),
                    upper=c(0.1,0.1,0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[8] <- bias
table_9$gompertz_CP[8] <- coverage





### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(3455354)
true_theta <- 4.02
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04                   
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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


table_9$exponential_MAE[9] <- bias
table_9$exponential_CP[9] <- coverage





### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(35567455)
true_theta <- 4.02
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04                   
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1, 0.5,1,0.5,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[9] <- bias
table_9$weibull_CP[9] <- coverage





### Frank copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(3532275)
true_theta <- 4.02
true_rho <- rho(frankCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04                   
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  frank_optim <- optim(c(0.1,0.1,0.1,0.1,2), frank_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,0.01),
                    upper=c(0.1,0.1,0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  frank_optim$par
  
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


table_9$gompertz_MAE[9] <- bias
table_9$gompertz_CP[9] <- coverage





### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(24521245)
true_theta <- 1.16
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.025                   #kidney failure
true_lambda2 <- 0.037                   #death

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

bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


table_9$exponential_MAE[10] <- bias
table_9$exponential_CP[10] <- coverage




### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(24008655)
true_theta <- 1.16
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.025                   #kidney failure
true_lambda2 <- 0.037                   #death

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
  gumbel_optim <- optim(c(0.2,0.2,0.2,0.2,2), gumbel_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,1.01),
                    upper=c(1,1,1,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[10] <- bias
table_9$weibull_CP[10] <- coverage



### Gumbel copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(982586609)
true_theta <- 1.16
true_rho <- rho(gumbelCopula(true_theta))
true_lambda1 <- 0.025                   #kidney failure
true_lambda2 <- 0.037                   #death

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
                    upper=c(0.2,0.2,0.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[10] <- bias
table_9$gompertz_CP[10] <- coverage





### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(24009545)
true_theta <- 1.26
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.60                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.84                   #death
true_b2 <- 0.04

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


table_9$exponential_MAE[11] <- bias
table_9$exponential_CP[11] <- coverage




### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(1211095)
true_theta <- 1.26
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.60                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.84                   #death
true_b2 <- 0.04

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
  gumbel_optim <- optim(c(0.2,0.2,0.2,0.2,2), gumbel_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,1.01),
                    upper=c(1,1,1,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[11] <- bias
table_9$weibull_CP[11] <- coverage





### Gumbel copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(985006545)
true_theta <- 1.26
true_rho <- rho(gumbelCopula(true_theta))
true_a1 <- 0.60                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.84                   #death
true_b2 <- 0.04

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
  
  g1_lw <- -0.1
  g1_up <- 1
  l1_lw <- 0.01
  l1_up <- 1
  g2_lw <- -0.1
  g2_up <- 1
  l2_lw <- 0.01
  l2_up <- 1
  t_lw <- 1.01
  t_up <- 10
  
  gumbel_optim <- optim(c(0.01, 0.04, 0.01, 0.03,1.17), gumbel_gompertz_loglik_sim, method="L-BFGS-B",lower=c(g1_lw, l1_lw, g2_lw, l2_lw, t_lw),
                    upper=c(g1_up, l1_up, g2_up, l2_up, t_up), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[11] <- bias
table_9$gompertz_CP[11] <- coverage





### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(240090895)
true_theta <- 1.17
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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


table_9$exponential_MAE[12] <- bias
table_9$exponential_CP[12] <- coverage





### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(8009395)
true_theta <- 1.17
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  gumbel_optim <- optim(c(0.2,0.2,0.2,0.2,2), gumbel_weibull_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,1.01),
                    upper=c(1,1,1,1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$weibull_MAE[12] <- bias
table_9$weibull_CP[12] <- coverage





### Gumbel copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(709909395)
true_theta <- 1.17
true_rho <- rho(gumbelCopula(true_theta))
true_gamma1 <- 0.001                   #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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
  gumbel_optim <- optim(c(0.01,0.04,0.01,0.03,2), gumbel_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.1,0.01,-0.1,0.01,1.01),
                    upper=c(0.2,0.2,0.2,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


table_9$gompertz_MAE[12] <- bias
table_9$gompertz_CP[12] <- coverage






### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(24070085)
true_theta <- 0.37
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.026                   #kidney failure
true_lambda2 <- 0.028                   #death

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


table_9$exponential_MAE[1] <- bias
table_9$exponential_CP[1] <- coverage




### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(240766332)
true_theta <- 0.37
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.026                   #kidney failure
true_lambda2 <- 0.028                   #death

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


table_9$weibull_MAE[1] <- bias
table_9$weibull_CP[1] <- coverage



### Normal copula ###
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(62213332)
true_theta <- 0.37
true_rho <- rho(normalCopula(true_theta))
true_lambda1 <- 0.026                   #kidney failure
true_lambda2 <- 0.028                   #death

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
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


table_9$gompertz_MAE[1] <- bias
table_9$gompertz_CP[1] <- coverage





### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Exponential
set.seed(244423085)
true_theta <- 0.46
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 0.61                    #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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


table_9$exponential_MAE[2] <- bias
table_9$exponential_CP[2] <- coverage




### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Weibull
set.seed(800966332)
true_theta <- 0.46
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 0.61                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.1),upper=c(1.5,0.2,2,0.2,0.9), 
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


table_9$weibull_MAE[2] <- bias
table_9$weibull_CP[2] <- coverage





### Normal copula ###
#True survival distribution = Weibull
#Assumed survival distribution = Gompertz
set.seed(7799332)
true_theta <- 0.46
true_rho <- rho(normalCopula(true_theta))
true_a1 <- 0.61                   #kidney failure
true_b1 <- 0.06
true_a2 <- 0.86                   #death
true_b2 <- 0.04

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


table_9$gompertz_MAE[2] <- bias
table_9$gompertz_CP[2] <- coverage





### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Exponential
set.seed(71121885)
true_theta <- 0.40
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.001                    #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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


table_9$exponential_MAE[3] <- bias
table_9$exponential_CP[3] <- coverage





### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Weibull
set.seed(5112132)
true_theta <- 0.40
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.001                    #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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


table_9$weibull_MAE[3] <- bias
table_9$weibull_CP[3] <- coverage





### Normal copula ###
#True survival distribution = Gompertz
#Assumed survival distribution = Gompertz
set.seed(76612132)
true_theta <- 0.40
true_rho <- rho(normalCopula(true_theta))
true_gamma1 <- 0.001                    #kidney failure
true_lambda1 <- 0.04
true_gamma2 <- 0.001                   #death
true_lambda2 <- 0.03

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

table_9$gompertz_MAE[3] <- bias
table_9$gompertz_CP[3] <- coverage


