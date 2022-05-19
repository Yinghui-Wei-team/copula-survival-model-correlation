#######################################################################################
#
#   Filename    :	      Functions.R    												  
#                                                                                     
#   Project     :       Article "Estimating the correlation between semi-competing risk survival endpoints"                                                             
#   Authors     :       L Sorrell, Y Wei, M Wojtys and P Rowe                                                               
#   Date        :       01/06/2021
#																				  
#   R Version   :       R-3.6.1                                                              
#
#   Required R packages :  copula, mvtnorm
#
#
########################################################################################
library(copula)
library(mvtnorm)

#Clayton loglikelihood with Exponential survival distributions
clayton_loglik <- function(para, X, Y, d1, d2){
  l1 <- para[1] #hazard rate for non-terminal event
  l2 <- para[2] #hazard rate for terminal event
  theta <- para[3] #association parameter
  
  C <- (exp(-l1*X)^(-theta)+exp(-l2*Y)^(-theta)-1)^(-1/theta) #copula
  
  part1 <- ifelse(d1*d2==1,(log(1+theta)+(1+2*theta)*log(C)-(theta+1)*log(exp(-l1*X))-(theta+1)*log(exp(-l2*Y))+log(l1)-l1*X+log(l2)-l2*Y),0) #both events
  part2 <- ifelse(d1*(1-d2)==1,((theta+1)*log(C)-(theta+1)*log(exp(-l1*X))+log(l1)-l1*X),0) #non-terminal event only
  part3 <- ifelse(((1-d1)*(d2))==1,((theta+1)*log(C)-(theta+1)*log(exp(-l2*Y))+log(l2)-l2*Y),0) #terminal event only
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0) #both events censored
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

#Frank loglikelihood with Exponential survival distributions
frank_loglik <- function(para, X, Y, d1, d2){
  l1 <- para[1]
  l2 <- para[2]
  theta <- para[3]
  
  C <- -1/theta * log(((1-exp(-theta)-(1-exp(-theta*exp(-l1*X)))*(1-exp(-theta*exp(-l2*Y)))))/(1-exp(-theta)))
  
  part1 <- ifelse(d1*d2==1,(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*exp(-l1*X))-1)-log(exp(theta*exp(-l2*Y))-1)+log(l1)-l1*X+log(l2)-l2*Y),0)
  part2 <- ifelse(d1*(1-d2)==1,(log((1-exp(theta*C))/(1-exp(theta*exp(-l1*X))))+log(l1)-l1*X),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log((1-exp(theta*C))/(1-exp(theta*exp(-l2*Y))))+log(l2)-l2*Y),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

#Gumbel loglikelihood with Exponential survival distributions
gumbel_loglik <- function(para, X, Y, d1, d2){
  l1 <- para[1]
  l2 <- para[2]
  theta <- para[3]
  
  C <- exp(-((-log(exp(-l1*X)))^(theta)+(-log(exp(-l2*Y)))^(theta))^(1/theta))
  
  part1 <- ifelse(d1*d2==1,(log(C)+(theta-1)*log(-log(exp(-l1*X)))+(theta-1)*log(-log(exp(-l2*Y)))+log(theta-1+((-log(exp(-l1*X)))^theta+(-log(exp(-l2*Y)))^theta)^(1/theta))-log(exp(-l1*X))-log(exp(-l2*Y))-(2*theta-1)*log(-log(C))+log(l1)-l1*X+log(l2)-l2*Y),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(C)+(theta-1)*log(-log(exp(-l1*X)))-log(exp(-l1*X))-(theta-1)*log(-log(C))+log(l1)-l1*X),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log(C)+(theta-1)*log(-log(exp(-l2*Y)))-log(exp(-l2*Y))-(theta-1)*log(-log(C))+log(l2)-l2*Y),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

#Normal loglikelihood with Exponential survival distributions
normal_loglik <- function(para, X, Y, d1, d2){
  l1 <- para[1]
  l2 <- para[2]
  rho <- para[3]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(pexp(X.1,l1) )*qnorm(pexp(Y.1,l2))-
                                                           rho^2*(qnorm(pexp(X.1,l1) )^2 + qnorm(pexp(Y.1,l2))^2)))/((2*(1-rho^2))))+
                                      log(l1)-l1*X.1 +log(l2)-l2*Y.1),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(pexp(Y.2,l2)), mean=rho*qnorm(pexp(X.2,l1)),
                                              sd=sqrt(1-rho^2), lower.tail=F)*l1*exp(-l1*X.2))),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(pexp(X.3,l1)), mean=rho*qnorm(pexp(Y.3,l2)),
                                              sd=sqrt(1-rho^2), lower.tail=F)*l2*exp(-l2*Y.3))),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(pexp(X.4,l1),pexp(Y.4,l2) )),1,normal_cdf,cov_matrix))),0)
  
  loglik <- (part1+part2+part3+part4) 
  return(loglik)
}


#Clayton loglikelihood with Gompertz survival distributions
clayton_gompertz_loglik <- function(para, X, Y, d1, d2){
  g1 <- para[1] #shape for non-terminal event
  l1 <- para[2] #rate for non-terminal event
  g2 <- para[3] #shape for terminal event
  l2 <- para[4] #rate for terminal event
  theta <- para[5] #association parameter
  
  S1 <- exp(-l1/g1*(exp(g1*X)-1))
  S2 <- exp(-l2/g2*(exp(g2*Y)-1))
  f1 <- l1*exp(g1*X-l1/g1*(exp(g1*X)-1))
  f2 <- l2*exp(g2*Y-l2/g2*(exp(g2*Y)-1))
  S1[which(S1 < 0.1^8)]<-0.1^8
  S2[which(S2 < 0.1^8)]<- 0.1^8
  f1[which(f1 < 0.1^8)]<-0.1^8
  f2[which(f2 < 0.1^8)]<- 0.1^8
  C <- (S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  
  part1 <- ifelse(d1*d2==1,(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(C^(1+theta)*f1)-log(S1^(1+theta))),0)
  part3 <- ifelse(d2*(1-d1)==1,(log(C^(1+theta)*f2)-log(S2^(1+theta))),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik<-sum(part1+part2+part3+part4) 
  return(loglik)
}


#Frank loglikelihood with Gompertz survival distributions
frank_gompertz_loglik <- function(para, X, Y, d1, d2){
  g1 <- para[1]
  l1 <- para[2]
  g2 <- para[3]
  l2 <- para[4]
  theta <- para[5]
  
  C <- -1/theta * log(((1-exp(-theta)-(1-exp(-theta*exp(-l1/g1*(exp(g1*X)-1))))*(1-exp(-theta*exp(-l2/g2*(exp(g2*Y)-1))))))/(1-exp(-theta)))
  
  part1 <- ifelse(d1*d2==1,(log(theta*exp(theta*C)*(exp(theta*C)-1)*l1*exp(g1*X-l1/g1*(exp(g1*X)-1))*l2*exp(g2*Y-l2/g2*(exp(g2*Y)-1)))-log((exp(theta*exp(-l1/g1*(exp(g1*X)-1)))-1)*(exp(theta*exp(-l2/g2*(exp(g2*Y)-1)))-1))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(((1-exp(theta*C))*l1*exp(g1*X-l1/g1*(exp(g1*X)-1)))/(1-exp(theta*exp(-l1/g1*(exp(g1*X)-1)))))),0) 
  part3 <- ifelse(d2*(1-d1)==1,(log(((1-exp(theta*C))*l2*exp(g2*Y-l2/g2*(exp(g2*Y)-1)))/(1-exp(theta*exp(-l2/g2*(exp(g2*Y)-1)))))),0) 
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}

#Gumbel loglikelihood with Gompertz survival distributions
gumbel_gompertz_loglik <- function(para, X, Y, d1, d2){
  g1 <- para[1]
  l1 <- para[2]
  g2 <- para[3]
  l2 <- para[4]
  theta <- para[5]
  
  C <- exp(-((-log(exp(-l1/g1*(exp(g1*X)-1))))^(theta)+(-log(exp(-l2/g2*(exp(g2*Y)-1))))^(theta))^(1/theta))
  
  part1 <- ifelse(d1*d2==1,(log(C*(-log(exp(-l1/g1*(exp(g1*X)-1))))^(theta-1)*(-log(exp(-l2/g2*(exp(g2*Y)-1))))^(theta-1)*(theta-1-log(C))*l1*exp(g1*X-l1/g1*(exp(g1*X)-1))*l2*exp(g2*Y-l2/g2*(exp(g2*Y)-1)))-log(exp(-l1/g1*(exp(g1*X)-1))*exp(-l2/g2*(exp(g2*Y)-1))*(-log(C))^(2*theta-1))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(C*(-log(exp(-l1/g1*(exp(g1*X)-1))))^(theta-1)*l1*exp(g1*X-l1/g1*(exp(g1*X)-1)))-log(exp(-l1/g1*(exp(g1*X)-1))*(-log(C))^(theta-1))),0)
  part3 <- ifelse(d2*(1-d1)==1,(log(C*(-log(exp(-l2/g2*(exp(g2*Y)-1))))^(theta-1)*l2*exp(g2*Y-l2/g2*(exp(g2*Y)-1)))-log(exp(-l2/g2*(exp(g2*Y)-1))*(-log(C))^(theta-1))),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)  
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}


#Normal loglikelihood with Gompertz survival distributions
normal_gompertz_loglik <- function(para, X, Y, d1, d2){
  g1 <- para[1]
  l1 <- para[2]
  g2 <- para[3]
  l2 <- para[4]
  rho <- para[5]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2    #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(1-exp(-l1/g1*(exp(g1*X.1)-1)))*qnorm(1-exp(-l2/g2*(exp(g2*Y.1)-1)))-
                                                           rho^2*(qnorm(1-exp(-l1/g1*(exp(g1*X.1)-1)))^2 + qnorm(1-exp(-l2/g2*(exp(g2*Y.1)-1)))^2)))/((2*(1-rho^2))))+
                                      log(l1*exp(g1*X.1-l1/g1*(exp(g1*X.1)-1)))+log(l2*exp(g2*Y.1-l2/g2*(exp(g2*Y.1)-1)))),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(1-exp(-l2/g2*(exp(g2*Y.2)-1))), mean=rho*qnorm(1-exp(-l1/g1*(exp(g1*X.2)-1))),
                                              sd=sqrt(1-rho^2), lower.tail=F)*l1*exp(g1*X.2-l1/g1*(exp(g1*X.2)-1)))),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(1-exp(-l1/g1*(exp(g1*X.3)-1))), mean=rho*qnorm(1-exp(-l2/g2*(exp(g2*Y.3)-1))),
                                              sd=sqrt(1-rho^2), lower.tail=F)*l2*exp(g2*Y.3-l2/g2*(exp(g2*Y.3)-1)))),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(1-exp(-l1/g1*(exp(g1*X.4)-1)),
                                                          1-exp(-l2/g2*(exp(g2*Y.4)-1)))),1,normal_cdf,cov_matrix))),0)
  
  
  loglik <- (part1+part2+part3+part4)
  return(loglik)
}


#Clayton loglikelihood with Weibull survival distributions
clayton_weibull_loglik <- function(para, X, Y, d1, d2){
  a1 <- para[1] #shape for non-terminal event
  b1 <- para[2] #scale for non-terminal event
  a2 <- para[3] #shape for terminal event
  b2 <- para[4] #scale for terminal event
  theta <- para[5] #association parameter
  
  S1 <- exp(-b1*X^a1)
  S2 <- exp(-b2*Y^a2)
  f1 <- b1*a1*X^(a1-1)*exp(-b1*X^a1) 
  f2 <- b2*a2*Y^(a2-1)*exp(-b2*Y^a2) 
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  
  C <- (S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  
  part1 <- ifelse(d1*d2==1,(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(C^(1+theta)*f1)-log(S1^(1+theta))),0)
  part3 <- ifelse(d2*(1-d1)==1,(log(C^(1+theta)*f2)-log(S2^(1+theta))),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik<-sum(part1+part2+part3+part4) 
  return(loglik)
}

#Frank loglikelihood with Weibull survival distributions
frank_weibull_loglik <- function(para, X, Y, d1, d2){
  a1 <- para[1]
  b1 <- para[2]
  a2 <- para[3]
  b2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-b1*X^a1)
  S2 <- exp(-b2*Y^a2)
  f1 <- b1*a1*X^(a1-1)*exp(-b1*X^a1) 
  f2 <- b2*a2*Y^(a2-1)*exp(-b2*Y^a2) 
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  
  C <- -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  
  part1 <- ifelse(d1*d2==1,(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))),0)
  part3 <- ifelse(d2*(1-d1)==1,(log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}

#Gumbel loglikelihood with Weibull survival distributions
gumbel_weibull_loglik <- function(para, X, Y, d1, d2){
  a1 <- para[1]
  b1 <- para[2]
  a2 <- para[3]
  b2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-b1*X^a1)
  S2 <- exp(-b2*Y^a2)
  f1 <- b1*a1*X^(a1-1)*exp(-b1*X^a1) 
  f2 <- b2*a2*Y^(a2-1)*exp(-b2*Y^a2) 
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  
  C <- exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  
  part1 <- ifelse(d1*d2==1,(log(C*(-log(S1))^(theta-1)*(-log(S2))^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1))),0)
  part2 <- ifelse(d1*(1-d2)==1,(log(C*(-log(S1))^(theta-1)*f1)-log(S1*(-log(C))^(theta-1))),0)
  part3 <- ifelse(d2*(1-d1)==1,(log(C*(-log(S2))^(theta-1)*f2)-log(S2*(-log(C))^(theta-1))),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}


#Normal loglikelihood with Weibull survival distributions
normal_weibull_loglik <- function(para, X, Y, d1, d2){
  a1 <- para[1]
  b1 <- para[2]
  a2 <- para[3]
  b2 <- para[4]
  rho <- para[5]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2    #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(1-exp(-b1*X.1^a1))*qnorm(1-exp(-b2*Y.1^a2))-
                                                           rho^2*(qnorm(1-exp(-b1*X.1^a1))^2 + qnorm(1-exp(-b2*Y.1^a2))^2)))/((2*(1-rho^2))))+
                                      log(b1*a1*X.1^(a1-1)*exp(-b1*X.1^a1))+log(b2*a2*Y.1^(a2-1)*exp(-b2*Y.1^a2))),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(1-exp(-b2*Y.2^a2)), mean=rho*qnorm(1-exp(-b1*X.2^a1)),
                                              sd=sqrt(1-rho^2), lower.tail=F)*(b1*a1*X.2^(a1-1)*exp(-b1*X.2^a1)))),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(1-exp(-b1*X.3^a1)), mean=rho*qnorm(1-exp(-b2*Y.3^a2)),
                                              sd=sqrt(1-rho^2), lower.tail=F)*(b2*a2*Y.3^(a2-1)*exp(-b2*Y.3^a2)))),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(1-exp(-b1*X.4^a1),1-exp(-b2*Y.4^a2))),1,normal_cdf,cov_matrix))),0)
  
  loglik <- (part1+part2+part3+part4)
  return(loglik);
}


##############################################################
#Repeat of functions above used in simulations with survival 
#functions and pdfs prevented from getting too small
##############################################################

clayton_loglik_sim <- function(para, X, Y, d1, d2){
  lambda1 <- para[1]
  lambda2 <- para[2]
  theta <- para[3]
  
  S1 <- exp(-lambda1*X)
  S2 <- exp(-lambda2*Y)
  f1 <- lambda1*exp(-lambda1*X)
  f2 <- lambda2*exp(-lambda2*Y)
  S1[which(S1 < 0.1^8)]<-0.1^8
  S2[which(S2 < 0.1^8)]<- 0.1^8
  f1[which(f1 < 0.1^8)]<-0.1^8
  f2[which(f2 < 0.1^8)]<- 0.1^8
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta)))
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  part4 <- ((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

clayton_weibull_loglik_sim <- function(para, X, Y, d1, d2){
  alpha1 <- para[1]
  beta1 <- para[2]
  alpha2 <- para[3]
  beta2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-beta1*X^alpha1)
  S2 <- exp(-beta2*Y^alpha2)
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta)))
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  part4 <- ((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

clayton_gompertz_loglik_sim <- function(para, X, Y, d1, d2){
  gamma1 <- para[1]
  lambda1 <- para[2]
  gamma2 <- para[3]
  lambda2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  S1[which(S1 < 0.1^8)]<-0.1^8
  S2[which(S2 < 0.1^8)]<- 0.1^8
  f1[which(f1 < 0.1^8)]<-0.1^8
  f2[which(f2 < 0.1^8)]<- 0.1^8
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta)))
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  part4 <- ((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

frank_loglik_sim <- function(para, X, Y, d1, d2){
  lambda1 <- para[1] #parameters to be estimated
  lambda2 <- para[2] 
  theta <- para[3]
  
  S1 <- exp(-lambda1*X) #pseudo X
  S2 <- exp(-lambda2*Y) #pseudo Y
  
  f1 <- lambda1*exp(-lambda1*X)
  f2 <- lambda2*exp(-lambda2*Y)
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
  part2 <- d1*(1-d2)*(log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))) 
  part3 <-((1-d1)*(d2))*(log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))) 
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

frank_weibull_loglik_sim <- function(para, X, Y, d1, d2){
  alpha1 <- para[1]
  beta1 <- para[2]
  alpha2 <- para[3]
  beta2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-beta1*X^alpha1)
  S2 <- exp(-beta2*Y^alpha2)
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  C[which(C<0.1^8)] <- 0.1^8
  part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
  part2 <- d1*(1-d2)*(log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))) 
  part3 <-((1-d1)*(d2))*(log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))) 
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

frank_gompertz_loglik_sim <- function(para, X, Y, d1, d2){
  gamma1 <- para[1]
  lambda1 <- para[2]
  gamma2 <- para[3]
  lambda2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  S1[which(S1 < 0.1^8)] <- 0.1^8
  S2[which(S2 < 0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  C[which(C<0.1^8)] <- 0.1^8
  part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
  part2 <- d1*(1-d2)*(log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))) 
  part3 <-((1-d1)*(d2))*(log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))) 
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

gumbel_loglik_sim <- function(para, X, Y, d1, d2){
  lambda1 <- para[1] 
  lambda2 <- para[2] 
  theta <- para[3]
  
  S1 <- exp(-lambda1*X)
  S2 <- exp(-lambda2*Y) 
  f1 <- lambda1*exp(-lambda1*X)
  f2 <- lambda2*exp(-lambda2*Y)
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C <- exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  C[which(C<0.1^8)] <- 0.1^8
  part1 <- d1*d2*(log(C*(-log(S1))^(theta-1)*(-log(S2))^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1)))
  part2 <- d1*(1-d2)*(log(C*(-log(S1))^(theta-1)*f1)-log(S1*(-log(C))^(theta-1)))
  part3 <-((1-d1)*(d2))*(log(C*(-log(S2))^(theta-1)*f2)-log(S2*(-log(C))^(theta-1)))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

gumbel_weibull_loglik_sim <- function(para, X, Y, d1, d2){
  alpha1 <- para[1]
  beta1 <- para[2]
  alpha2 <- para[3]
  beta2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-beta1*X^alpha1)
  S2 <- exp(-beta2*Y^alpha2)
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C <- exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  part1 <- d1*d2*(log(C*(-log(S1))^(theta-1)*(-log(S2))^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1)))
  part2 <- d1*(1-d2)*(log(C*(-log(S1))^(theta-1)*f1)-log(S1*(-log(C))^(theta-1)))
  part3 <-((1-d1)*(d2))*(log(C*(-log(S2))^(theta-1)*f2)-log(S2*(-log(C))^(theta-1)))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

gumbel_gompertz_loglik_sim <- function(para, X, Y, d1, d2){
  gamma1 <- para[1]
  lambda1 <- para[2]
  gamma2 <- para[3]
  lambda2 <- para[4]
  theta <- para[5]
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  C <- exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  part1 <- d1*d2*(log(C*(-log(S1))^(theta-1)*(-log(S2))^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1)))
  part2 <- d1*(1-d2)*(log(C*(-log(S1))^(theta-1)*f1)-log(S1*(-log(C))^(theta-1)))
  part3 <-((1-d1)*(d2))*(log(C*(-log(S2))^(theta-1)*f2)-log(S2*(-log(C))^(theta-1)))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

normal_loglik_sim <- function(para, X, Y, d1, d2){
  lambda1 <- para[1]
  lambda2 <- para[2]
  rho <- para[3]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  p1 <- exp(-lambda1*X.1)
  p2 <- exp(-lambda2*Y.1)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  S1.1 <- 1-p1
  S2.1 <- 1-p2
  S1.1[which(S1.1<0.1^8)] <- 0.1^8
  S2.1[which(S2.1<0.1^8)] <- 0.1^8
  f1.1 <- lambda1*exp(-lambda1*X.1) 
  f2.1 <- lambda2*exp(-lambda2*Y.1) 
  f1.1[which(f1.1<0.1^8)] <- 0.1^8
  f2.1[which(f2.1<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1*X.1)
  p2 <- exp(-lambda2*Y.1)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  S1.2 <- 1-p1
  S2.2 <- 1-p2
  S1.2[which(S1.1<0.1^8)] <- 0.1^8
  S2.2[which(S2.1<0.1^8)] <- 0.1^8
  f1.2 <- lambda1*exp(-lambda1*X.1) 
  f2.2 <- lambda2*exp(-lambda2*Y.1) 
  f1.2[which(f1.1<0.1^8)] <- 0.1^8
  f2.2[which(f2.1<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1*X.1)
  p2 <- exp(-lambda2*Y.1)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  S1.3 <- 1-p1
  S2.3 <- 1-p2
  S1.3[which(S1.1<0.1^8)] <- 0.1^8
  S2.3[which(S2.1<0.1^8)] <- 0.1^8
  f1.3 <- lambda1*exp(-lambda1*X.1) 
  f2.3 <- lambda2*exp(-lambda2*Y.1) 
  f1.3[which(f1.1<0.1^8)] <- 0.1^8
  f2.3[which(f2.1<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1*X.4)
  p2 <- exp(-lambda2*Y.4)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  S1.4 <- 1-p1
  S2.4 <- 1-p2
  S1.4[which(S1.4<0.1^8)] <- 0.1^8
  S2.4[which(S2.4<0.1^8)] <- 0.1^8
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(S1.1)*qnorm(S2.2)-
                                                           rho^2*(qnorm(S1.1)^2 + qnorm(S2.2)^2)))/((2*(1-rho^2))))+
                                      log(f1.1)+log(f2.2)),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(S2.2), mean=rho*qnorm(S1.1),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f1.2)),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(S1.3), mean=rho*qnorm(S2.3),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f2.3)),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(S1.4,S2.4)),1,normal_cdf,cov_matrix))),0)
  
  loglik <- (part1+part2+part3+part4) 
  return(loglik)
}

normal_weibull_loglik_sim <- function(para, X, Y, d1, d2){
  alpha1 <- para[1]
  beta1 <- para[2]
  alpha2 <- para[3]
  beta2 <- para[4]
  rho <- para[5]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2    #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  p1 <- exp(-beta1*X.1^alpha1)
  p2 <- exp(-beta2*Y.1^alpha2)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.1 <- 1-p1
  X2.1 <- 1-p2
  f1.1 <- beta1*alpha1*X.1^(alpha1-1)*exp(-beta1*X.1^alpha1) 
  f2.1 <- beta2*alpha2*Y.1^(alpha2-1)*exp(-beta2*Y.1^alpha2) 
  f1.1[which(f1.1<0.1^8)] <- 0.1^8
  f2.1[which(f2.1<0.1^8)] <- 0.1^8
  
  p1 <- exp(-beta1*X.2^alpha1)
  p2 <- exp(-beta2*Y.2^alpha2)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.2 <- 1-p1
  X2.2 <- 1-p2
  f1.2 <- beta1*alpha1*X.2^(alpha1-1)*exp(-beta1*X.2^alpha1) 
  f2.2 <- beta2*alpha2*Y.2^(alpha2-1)*exp(-beta2*Y.2^alpha2) 
  f1.2[which(f1.2<0.1^8)] <- 0.1^8
  f2.2[which(f2.2<0.1^8)] <- 0.1^8
  
  p1 <- exp(-beta1*X.3^alpha1)
  p2 <- exp(-beta2*Y.3^alpha2)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.3 <- 1-p1
  X2.3 <- 1-p2
  f1.3 <- beta1*alpha1*X.3^(alpha1-1)*exp(-beta1*X.3^alpha1)
  f2.3 <- beta2*alpha2*Y.3^(alpha2-1)*exp(-beta2*Y.3^alpha2) 
  f1.3[which(f1.3<0.1^8)] <- 0.1^8
  f2.3[which(f2.3<0.1^8)] <- 0.1^8
  
  p1 <- exp(-beta1*X.4^alpha1)
  p2 <- exp(-beta2*Y.4^alpha2)
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.4 <- 1-p1
  X2.4 <- 1-p2
  
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(X1.1)*qnorm(X2.1)-
                                                           rho^2*(qnorm(X1.1)^2 + qnorm(X2.1)^2)))/((2*(1-rho^2))))+
                                      log(f1.1)+log(f2.1)),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(X2.2), mean=rho*qnorm(X1.2),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f1.2)),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(X1.3), mean=rho*qnorm(X2.3),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f2.3)),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(X1.4,X2.4)),1,normal_cdf,cov_matrix))),0)
  
  loglik <- (part1+part2+part3+part4)
  return(loglik);
}

normal_gompertz_loglik_sim <- function(para, X, Y, d1, d2){
  gamma1 <- para[1]
  lambda1 <- para[2]
  gamma2 <- para[3]
  lambda2 <- para[4]
  rho <- para[5]
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2    #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  p1 <- exp(-lambda1/gamma1*(exp(gamma1*X.1)-1))
  p2 <- exp(-lambda2/gamma2*(exp(gamma2*Y.1)-1))
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.1 <- 1-p1
  X2.1 <- 1-p2
  f1.1 <- lambda1*exp(gamma1*X.1-lambda1/gamma1*(exp(gamma1*X.1)-1))
  f2.1 <- lambda2*exp(gamma2*Y.1-lambda2/gamma2*(exp(gamma2*Y.1)-1))
  f1.1[which(f1.1<0.1^8)] <- 0.1^8
  f2.1[which(f2.1<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1/gamma1*(exp(gamma1*X.2)-1))
  p2 <- exp(-lambda2/gamma2*(exp(gamma2*Y.2)-1))
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.2 <- 1-p1
  X2.2 <- 1-p2
  f1.2 <- lambda1*exp(gamma1*X.2-lambda1/gamma1*(exp(gamma1*X.2)-1))
  f2.2 <- lambda2*exp(gamma2*Y.2-lambda2/gamma2*(exp(gamma2*Y.2)-1))
  f1.2[which(f1.2<0.1^8)] <- 0.1^8
  f2.2[which(f2.2<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1/gamma1*(exp(gamma1*X.3)-1))
  p2 <- exp(-lambda2/gamma2*(exp(gamma2*Y.3)-1))
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.3 <- 1-p1
  X2.3 <- 1-p2
  f1.3 <- lambda1*exp(gamma1*X.3-lambda1/gamma1*(exp(gamma1*X.3)-1))
  f2.3 <- lambda2*exp(gamma2*Y.3-lambda2/gamma2*(exp(gamma2*Y.3)-1))
  f1.3[which(f1.3<0.1^8)] <- 0.1^8
  f2.3[which(f2.3<0.1^8)] <- 0.1^8
  
  p1 <- exp(-lambda1/gamma1*(exp(gamma1*X.4)-1))
  p2 <- exp(-lambda2/gamma2*(exp(gamma2*Y.4)-1))
  p1[which(p1<0.1^8)] <- 0.1^8
  p2[which(p2<0.1^8)] <- 0.1^8
  X1.4 <- 1-p1
  X2.4 <- 1-p2
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(X1.1)*qnorm(X2.1)-
                                                           rho^2*(qnorm(X1.1)^2 + qnorm(X2.1)^2)))/((2*(1-rho^2))))+
                                      log(f1.1)+log(f2.1)),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(X2.2), mean=rho*qnorm(X1.2),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f1.2)),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(X1.3), mean=rho*qnorm(X2.3),
                                              sd=sqrt(1-rho^2), lower.tail=F)*f2.3)),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(X1.4,X2.4)),1,normal_cdf,cov_matrix))),0)
  
  
  loglik <- (part1+part2+part3+part4)
  return(loglik)
}
