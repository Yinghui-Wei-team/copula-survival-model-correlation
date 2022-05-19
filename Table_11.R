#######################################################################################
#
#   Filename    :	      Table 11.R    												  
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

#Prepare table
exponential_MAE <- rep(NA, 20)
exponential_CP <- rep(NA, 20)
exponential_boundary <- rep(NA, 20)
weibull_MAE <- rep(NA, 20)
weibull_CP <- rep(NA, 20)
weibull_boundary <- rep(NA, 20)
gompertz_MAE <- rep(NA, 20)
gompertz_CP <- rep(NA, 20)
gompertz_boundary <- rep(NA, 20)
rho <- c(rep(0,4), rep(0.25,4), rep(0.5, 4), rep(0.75,4), rep(0.9,4))
copula <- rep(c("Normal", "Clayton", "Frank", "Gumbel"),5)
table_11 <- data.frame(rho, copula, exponential_MAE, exponential_CP, exponential_boundary,
                       weibull_MAE, weibull_CP,  weibull_boundary,
                       gompertz_MAE, gompertz_CP,  gompertz_boundary)

true_lambda1 <- log(2)/5           #hazard rate for time-to-progression
true_lambda2 <- log(2)/11          #hazard rate for death
n <- 100                           #number of individuals 
runs <- 1000                       #number of data sets

#Rho=0, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(421439800)
true_theta <- 0.001
true_rho <- rho(claytonCopula(true_theta)) 

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)

for(i in 1:(runs)){
  
  for(k in 1:n){   #loop to simulate U an V
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
  clayton_optim <- optim(c(0.1,0.1,1), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,12.2), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-clayton_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- clayton_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- clayton_optim$par[3]-1.96*se[3] #lower ci for theta
  
  if(lower_ci_theta < -1) {lower_ci_theta <- -1}
  if(lower_ci_theta=="NaN") {next}  

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


#Report
table_11$exponential_MAE[2] <- bias
table_11$exponential_CP[2] <- coverage




#Rho=0, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(653266467)
true_theta <- 0.001
true_rho <- rho(claytonCopula(true_theta)) 

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
  if(lower_ci_theta < -1) {lower_ci_theta <- -1}
  
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


#Report
table_11$weibull_MAE[2] <- bias
table_11$weibull_CP[2] <- coverage



#Rho=0, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(866775645)
true_theta <- 0.001
true_rho <- rho(claytonCopula(true_theta)) 

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
  if(lower_ci_theta < -1) {lower_ci_theta <- -1}
  if(lower_ci_theta=="NaN") {next}  
  
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


#Report
table_11$gompertz_MAE[2] <- bias
table_11$gompertz_CP[2] <- coverage




#Rho=0, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(908994466)
true_theta <- 0.001
true_rho <- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.1,0.1,1), frank_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),upper=c(0.3,0.3,23.9), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
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


#Report
table_11$exponential_MAE[3] <- bias
table_11$exponential_CP[3] <- coverage



#Rho=0, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(216664553)
true_theta <- 0.01
true_rho <- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.7,0.23,0.7,0.1,1), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.4,0.5,1.4,0.2,5), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$weibull_MAE[3] <- bias
table_11$weibull_CP[3] <- coverage




#Rho=0, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(36696545)
true_theta <- 5.8
true_rho <- rho(frankCopula(true_theta))

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
                    upper=c(0.1,0.3,0.1,0.1,7), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$gompertz_MAE[3] <- bias
table_11$gompertz_CP[3] <- coverage





#Rho=0, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(231351212)
true_theta <- 1.01
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
counterless=0
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
  gumbel_optim <- optim(c(0.02,0.02,1.05), gumbel_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,1.01),upper=c(0.2,0.1,7),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if (lower_ci_theta <1) {lower_ci_theta=1}
  if (lower_ci_theta==1) {counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$exponential_MAE[4] <- bias
table_11$exponential_CP[4] <- coverage
table_11$exponential_boundary[4] <- boundary



#Rho=0, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(5455409)
true_theta <- 1.01
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
counterless=0
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
  if(lower_ci_theta==1) {counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$weibull_MAE[4] <- bias
table_11$weibull_CP[4] <- coverage
table_11$weibull_boundary[4] <- boundary




#Rho=0, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(982586851)
true_theta <- 1.01
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
counterless=0
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
                    upper=c(0.2,0.3,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  if(lower_ci_theta==1){counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$gompertz_MAE[4] <- bias
table_11$gompertz_CP[4] <- coverage
table_11$gompertz_boundary[4] <- boundary




#Rho=0, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(290406765)
true_theta <- 0
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.4), normal_loglik, method="L-BFGS-B",lower=c(0.01,0.01,0.01),upper=c(0.4,0.4,0.97), 
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


#Report
table_11$exponential_MAE[1] <- bias
table_11$exponential_CP[1] <- coverage




#Rho=0, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(814446646)
true_theta <- 0
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.001),upper=c(2,0.3,3,0.2,0.95), 
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


#Report
table_11$weibull_MAE[1] <- bias
table_11$weibull_CP[1] <- coverage



#Rho=0, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(902133332)
true_theta <- 0.91
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.03,0.01,-0.05,0.01,0.001),upper=c(0.2,0.3,0.9,0.3,0.97), 
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


#Report
table_11$gompertz_MAE[1] <- bias
table_11$gompertz_CP[1] <- coverage




#Rho=0.25, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(42198009)
true_theta <- 0.4
true_rho <- rho(claytonCopula(true_theta)) 

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
  clayton_optim <- optim(c(0.3,0.3,2), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),
                    upper=c(0.3,0.3,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
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


#Report
table_11$exponential_MAE[6] <- bias
table_11$exponential_CP[6] <- coverage



#Rho=0.25, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(653221167)
true_theta <- 0.4
true_rho <- rho(claytonCopula(true_theta)) 

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


#Report
table_11$weibull_MAE[6] <- bias
table_11$weibull_CP[6] <- coverage





#Rho=0.25, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(86677457)
true_theta <- 0.4
true_rho <- rho(claytonCopula(true_theta)) 

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


#Report
table_11$gompertz_MAE[6] <- bias
table_11$gompertz_CP[6] <- coverage






#Rho=0.25, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(90899066)
true_theta <- 1.55
true_rho <- rho(frankCopula(true_theta))

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


#Report
table_11$exponential_MAE[7] <- bias
table_11$exponential_CP[7] <- coverage




#Rho=0.25, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(212321957)
true_theta <- 1.55
true_rho <- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.5,2.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$weibull_MAE[7] <- bias
table_11$weibull_CP[7] <- coverage




#Rho=0.25, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(36696545)
true_theta <- 5.8
true_rho <- rho(frankCopula(true_theta))

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
                    upper=c(0.1,0.3,0.1,0.1,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$gompertz_MAE[7] <- bias
table_11$gompertz_CP[7] <- coverage




#Rho=0.25, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(231351212)
true_theta <- 1.20
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
counterless=0

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
  gumbel_optim <- optim(c(0.02,0.02,1.01), gumbel_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,1.01),upper=c(0.3,0.3,6.9), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #optimise
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if (lower_ci_theta <1) {lower_ci_theta=1}
  if (lower_ci_theta==1) {counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$exponential_MAE[8] <- bias
table_11$exponential_CP[8] <- coverage
table_11$exponential_boundary[8] <- boundary




#Rho=0.25, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(80554509)
true_theta <- 1.20
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
counterless=0
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
  if(lower_ci_theta==1){counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$weibull_MAE[8] <- bias
table_11$weibull_CP[8] <- coverage
table_11$weibull_boundary[8] <- boundary





#Rho=0.25, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(982586800)
true_theta <- 1.20
true_rho <- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
counterless=0

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
                    upper=c(0.2,0.3,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  if(lower_ci_theta==1){counterless=counterless+1}
  
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
boundary <- counterless

#Report
table_11$gompertz_MAE[8] <- bias
table_11$gompertz_CP[8] <- coverage
table_11$gompertz_boundary[8] <- boundary




#Rho=0.25, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(29042458)
true_theta <- 0.26
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.01,0.01,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(0.3,0.3,0.9),
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

#Report
table_11$exponential_MAE[5] <- bias
table_11$exponential_CP[5] <- coverage




#Rho=0.25, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(81211065)
true_theta <- 0.26
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),upper=c(2,0.3,3,0.2,0.9), 
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


#Report
table_11$weibull_MAE[5] <- bias
table_11$weibull_CP[5] <- coverage





#Rho=0.25, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(90098998)
true_theta <- 0.26
true_rho <- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.03,0.01,-0.05,0.01,0.01),upper=c(0.2,0.3,0.9,0.3,0.9), 
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


#Report
table_11$gompertz_MAE[5] <- bias
table_11$gompertz_CP[5] <- coverage




#Rho=0.5, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(421967429)
true_theta <- 1.06
true_rho<- rho(claytonCopula(true_theta)) #try 0.5

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
  clayton_optim <- optim(c(0.3,0.3,2), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),
                    upper=c(0.3,0.3,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
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
}

#Results 
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


#Report
table_11$exponential_MAE[10] <- bias
table_11$exponential_CP[10] <- coverage




#Rho=0.5, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(65300907)
true_theta <- 1.06
true_rho<- rho(claytonCopula(true_theta)) #try 0.5

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


#Report
table_11$weibull_MAE[10] <- bias
table_11$weibull_CP[10] <- coverage





#Rho=0.5, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(8660907)
true_theta <- 1.06
true_rho<- rho(claytonCopula(true_theta)) #try 0.5

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


#Report
table_11$gompertz_MAE[10] <- bias
table_11$gompertz_CP[10] <- coverage




#Rho=0.5, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(9033066)
true_theta <- 3.45
true_rho <- rho(frankCopula(true_theta))

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


#Report
table_11$exponential_MAE[11] <- bias
table_11$exponential_CP[11] <- coverage




#Rho=0.5, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(2180980)
true_theta <- 3.45
true_rho <- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.5,2.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$weibull_MAE[11] <- bias
table_11$weibull_CP[11] <- coverage





#Rho=0.5, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(366090545)
true_theta <- 3.45
true_rho<- rho(frankCopula(true_theta))

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

#Report
table_11$gompertz_MAE[11] <- bias
table_11$gompertz_CP[11] <- coverage




#Rho=0.5, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(77874232)
true_theta <- 1.54
true_rho<- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
counterless=0

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
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,control=list(fnscale=-1),hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[3]+1.96*se[3] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[3]-1.96*se[3] #lower ci for theta
  if (lower_ci_theta <1) {lower_ci_theta=1}
  if (lower_ci_theta==1) {counterless=counterless+1}
  
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


#Report
table_11$exponential_MAE[12] <- bias
table_11$exponential_CP[12] <- coverage
table_11$exponential_boundary[12] <- counterless



#Rho=0.5, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(805332509)
true_theta <- 1.54
true_rho<- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
counterless=0

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
  if(lower_ci_theta==1) {counterless=counterless+1}
  
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


#Report
table_11$weibull_MAE[12] <- bias
table_11$weibull_CP[12] <- coverage
table_11$weibull_boundary[12] <- counterless



#Rho=0.5, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(121135003)
true_theta <- 1.54
true_rho<- rho(gumbelCopula(true_theta))

estimated_rho <- rep(0, runs)
lower_ci <- rep(0,runs)
upper_ci <-rep(0,runs)
estimated_theta <- rep(0,runs)
true_rho_vec <- rep(true_rho,runs)
counter_rho=0
U1 <- rep(0,n)
V1 <- rep(0,n)
counterless=0

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
                    upper=c(0.2,0.3,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  if(lower_ci_theta==1) {counterless=counterless+1}
  
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


#Report
table_11$gompertz_MAE[12] <- bias
table_11$gompertz_CP[12] <- coverage
table_11$gompertz_boundary[12] <- counterless




#Rho=0.5, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(290424421)
true_theta <- 0.52
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.01,0.01,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(0.3,0.3,0.9),
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


#Report
table_11$exponential_MAE[9] <- bias
table_11$exponential_CP[9] <- coverage




#Rho=0.5, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(812142246)
true_theta <- 0.52
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),upper=c(2,0.3,3,0.2,0.9), 
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


#Report
table_11$weibull_MAE[9] <- bias
table_11$weibull_CP[9] <- coverage





#Rho=0.5, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(90212998)
true_theta <- 0.52
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.03,0.01,-0.05,0.01,0.01),upper=c(0.2,0.3,0.9,0.3,0.9), 
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


#Report
table_11$gompertz_MAE[9] <- bias
table_11$gompertz_CP[9] <- coverage





#Rho=0.75, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(421965429)
true_theta <- 2.58
true_rho<- rho(claytonCopula(true_theta)) 

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
  clayton_optim <- optim(c(0.3,0.3,2), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),
                         upper=c(0.3,0.3,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                         hessian=TRUE)
  
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


#Report
table_11$exponential_MAE[14] <- bias
table_11$exponential_CP[14] <- coverage




#Rho=0.75, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(656653907)
true_theta <- 2.58
true_rho<- rho(claytonCopula(true_theta)) #try 0.75

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


#Report
table_11$weibull_MAE[14] <- bias
table_11$weibull_CP[14] <- coverage





#Rho=0.75, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(866012190)
true_theta <- 2.58
true_rho<- rho(claytonCopula(true_theta)) 

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


#Report
table_11$gompertz_MAE[14] <- bias
table_11$gompertz_CP[14] <- coverage




#Rho=0.75, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(9034121)
true_theta <- 6.71
true_rho<- rho(frankCopula(true_theta))

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


#Report
table_11$exponential_MAE[15] <- bias
table_11$exponential_CP[15] <- coverage




#Rho=0.75, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(213180980)
true_theta <- 6.71
true_rho<- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.01,0.01),
                    upper=c(1.5,0.5,2.5,0.2,12), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$weibull_MAE[15] <- bias
table_11$weibull_CP[15] <- coverage




#Rho=0.75, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(36609005)
true_theta <- 6.71
true_rho<- rho(frankCopula(true_theta))

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
                    upper=c(0.2,0.3,0.1,0.3,15), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$gompertz_MAE[15] <- bias
table_11$gompertz_CP[15] <- coverage




#Rho=0.75, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(778090023)
true_theta <- 2.29
true_rho<- rho(gumbelCopula(true_theta))

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
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
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


#Report
table_11$exponential_MAE[16] <- bias
table_11$exponential_CP[16] <- coverage




#Rho=0.75, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(801213250)
true_theta <- 2.29
true_rho<- rho(gumbelCopula(true_theta))

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


#Report
table_11$weibull_MAE[16] <- bias
table_11$weibull_CP[16] <- coverage





#Rho=0.75, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(121095003)
true_theta <- 2.29
true_rho<- rho(gumbelCopula(true_theta))

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
                    upper=c(0.2,0.3,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$gompertz_MAE[16] <- bias
table_11$gompertz_CP[16] <- coverage




#Rho=0.75, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(290400908)
true_theta <- 0.77
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.01,0.01,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(0.3,0.3,0.99),
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


#Report
table_11$exponential_MAE[13] <- bias
table_11$exponential_CP[13] <- coverage




#Rho=0.75, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(812333465)
true_theta <- 0.77
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),upper=c(2,0.3,3,0.2,0.95), 
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


#Report
table_11$weibull_MAE[13] <- bias
table_11$weibull_CP[13] <- coverage





#Rho=0.75, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(902133398)
true_theta <- 0.77
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.03,0.01,-0.05,0.01,0.01),upper=c(0.2,0.3,0.9,0.3,0.9), 
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


#Report
table_11$gompertz_MAE[13] <- bias
table_11$gompertz_CP[13] <- coverage





#Rho=0.9, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(221196542)
true_theta <- 5.55
true_rho<- rho(claytonCopula(true_theta)) 

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
  clayton_optim <- optim(c(0.3,0.3,2), clayton_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,0.01),
                    upper=c(0.3,0.3,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
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
}

#Results 
bias <- mean(abs(estimated_rho-true_rho_vec))
mse <- mean((estimated_rho-true_rho_vec)^2)
coverage <- (counter_rho / runs) * 100


#Report
table_11$exponential_MAE[18] <- bias
table_11$exponential_CP[18] <- coverage




#Rho=0.9, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(1121653907)
true_theta <- 5.55
true_rho<- rho(claytonCopula(true_theta)) #try 0.9

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


#Report
table_11$weibull_MAE[18] <- bias
table_11$weibull_CP[18] <- coverage





#Rho=0.9, Clayton copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(80091907)
true_theta <- 5.55
true_rho<- rho(claytonCopula(true_theta)) #try 0.9

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


#Report
table_11$gompertz_MAE[18] <- bias
table_11$gompertz_CP[18] <- coverage





#Rho=0.9, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(96654121)
true_theta <- 12.25
true_rho<- rho(frankCopula(true_theta))

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


#Report
table_11$exponential_MAE[19] <- bias
table_11$exponential_CP[19] <- coverage




#Rho=0.9, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(2132111980)
true_theta <- 12.25
true_rho<- rho(frankCopula(true_theta))

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
  frank_optim <- optim(c(0.2,0.2,0.2,0.2,2), frank_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),
                    upper=c(1.5,0.5,2.5,0.2,20), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$weibull_MAE[19] <- bias
table_11$weibull_CP[19] <- coverage





#Rho=0.9, Frank copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(366044900)
true_theta <- 12.25
true_rho<- rho(frankCopula(true_theta))

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
                    upper=c(0.1,0.3,0.1,0.2,17), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
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


#Report
table_11$gompertz_MAE[19] <- bias
table_11$gompertz_CP[19] <- coverage




#Rho=0.9, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(775551213)
true_theta <- 3.74
true_rho<- rho(gumbelCopula(true_theta))

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
  gumbel_optim <- optim(c(0.1,0.1,1), gumbel_loglik_sim, method="L-BFGS-B",lower=c(0.001,0.001,1),upper=c(0.3,0.3,12),
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


#Report
table_11$exponential_MAE[20] <- bias
table_11$exponential_CP[20] <- coverage




#Rho=0.9, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(82221325)
true_theta <- 3.74
true_rho<- rho(gumbelCopula(true_theta))

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


#Report
table_11$weibull_MAE[20] <- bias
table_11$weibull_CP[20] <- coverage





#Rho=0.9, Gumbel copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(121090098)
true_theta <- 3.74
true_rho<- rho(gumbelCopula(true_theta))

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
                    upper=c(0.2,0.3,0.5,0.2,10), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, control=list(fnscale=-1), 
                    hessian=TRUE)
  
  #Confidence Intervals
  fisher_info <- solve(-gumbel_optim$hessian)   #Fisher's Information Matrix
  se <- sqrt(diag(fisher_info)) #standard error
  
  upper_ci_theta <- gumbel_optim$par[5]+1.96*se[5] #upper ci for theta
  lower_ci_theta <- gumbel_optim$par[5]-1.96*se[5] #lower ci for theta
  if(lower_ci_theta<1){lower_ci_theta=1}
  if(lower_ci_theta==1) {counterless=counterless+1}
  
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


#Report
table_11$gompertz_MAE[20] <- bias
table_11$gompertz_CP[20] <- coverage





#Rho=0.9, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Exponential
set.seed(29040558)
true_theta <- 0.91
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.01,0.01,0.1), normal_loglik, method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(0.3,0.3,0.99),
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


#Report
table_11$exponential_MAE[17] <- bias
table_11$exponential_CP[17] <- coverage




#Rho=0.9, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Weibull
set.seed(814443465)
true_theta <- 0.91
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.1,0.1,0.1,0.1,0.4), normal_weibull_loglik_sim, method="L-BFGS-B",lower=c(0.01,0.01,0.01,0.001,0.01),upper=c(2,0.3,3,0.2,0.95), 
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


#Report
table_11$weibull_MAE[17] <- bias
table_11$weibull_CP[17] <- coverage





#Rho=0.9, Normal copula
#True survival distribution = Exponential
#Assumed survival distribution = Gompertz
set.seed(902133332)
true_theta <- 0.91
true_rho<- rho(normalCopula(true_theta))

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
  normal_optim <- optim(c(0.05, 0.05, 0.23, 0.02, 0.59), normal_gompertz_loglik_sim, method="L-BFGS-B",lower=c(-0.03,0.01,-0.05,0.01,0.01),upper=c(0.2,0.3,0.9,0.3,0.97), 
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


#Report
table_11$gompertz_MAE[17] <- bias
table_11$gompertz_CP[17] <- coverage

