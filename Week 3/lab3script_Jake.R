library(tidyverse)
covid_districts<- read_csv(file = "https://people.bath.ac.uk/kai21/ASI/data/local_authority_covid_deaths.csv")
covid_districts %>% select(local_authority,death_rate) -> covid_districts #select only the area and number of deaths

#### Question 1 ####
N=10000 #number of samples required
samples<-matrix(0,nrow=N,ncol=311) #matrix to hold samples
for (i in 1:N){ #for loop to generate samples
  samples[i,]<-rnbinom(311,size=10,prob=0.05) #call rnbinom() to give us samples
} 

#### Question 2 - Correct! ####
loglikelihood_nbinom<-function(size,prob,data){
  
  sum_log_dens <- function(x,size,prob){
    sum(dnbinom(x,size,prob,log = TRUE))  
  }
  
  map2_dbl(.x = size,
           .y = prob,
           .f = sum_log_dens,
           x = data)
  
}
loglikelihood_nbinom3<- function(logsize,logitprob,data){
  
  
  size <- exp(logsize)
  prob  <- exp(logitprob)/(1+exp(logitprob))
  
  loglikelihood_nbinom(size,prob,data)
  
}
negloglik_fn<-function(theta = c(0,0),
                       data = 1){
  
  -loglikelihood_nbinom3(logsize = theta[1],
                         logitprob = theta[2],
                         data    = data)
}
#Custom log-likelihood function in this case, not given in notes
loglik_expr<-expression(lgamma(y+exp(theta1))-lgamma(exp(theta1))-lgamma(y+1)+exp(theta1)*(theta2-log(1+exp(theta2)))-y*log(1+exp(theta2)))
loglik_deriv <- deriv(expr         = loglik_expr,
                         namevec      = c("theta1","theta2"),
                         function.arg = c("theta1","theta2","y"),
                         hessian      = TRUE)
#New function to compute gradient of negative log-likelihood
negloglik_grad <- function(theta = c(0,0),
                           data = 1){
  
  aux  <- loglik_deriv(theta1 = theta[1],
                          theta2 = theta[2],
                          y      = data)
  
  grad <- apply(attr(aux,"gradient"),2,sum)
  
  -grad
}
#New function to compute Hessian
negloglik_hess <- function(theta = c(0,0),
                           data = 1){
  
  aux  <- loglik_deriv(theta1 = theta[1],
                          theta2 = theta[2],
                          y      = data)
  
  hess<-apply(attr(aux,"hessian"),c(2,3),sum)
  
  -hess
}


parameter_values<-matrix(0,nrow=N,ncol=2) #to store MLEs
theta_star<-c(log(10),log(0.05/0.95))
grad_values<-matrix(0,nrow=N,ncol=2) #to store gradient at theta_star
hessian_values<-array(0,dim=c(2,2,N)) #store Hessian at theta_star
newton_direction_values<-matrix(0,nrow=N,ncol=2) #store Newton descent directions at theta_star
for (i in 1:N){
  run<-optim(par=c(0,0),fn=negloglik_fn,gr=negloglik_grad,data=samples[i,],
             method="BFGS",hessian=TRUE) #perform optimisation
  parameter_values[i,]<-run$par #store MLEs
  grad_values[i,]<-negloglik_grad(theta_star,samples[i,]) #store gradient at theta_star
  hessian_values[,,i]<-negloglik_hess(theta_star,samples[i,]) #store hessian at theta_star
  newton_direction_values[i,]<- -solve(hessian_values[,,i])%*%grad_values[i,] #calculate newton direction
}

#Question 6 
#Parts 1 and 2
theta_star_matrix<-matrix(0,nrow=N,ncol=2)
theta_star_matrix[,1]<-rep(theta_star[1],N)
theta_star_matrix[,2]<-rep(theta_star[2],N)
differences<-parameter_values-theta_star_matrix
varcov_differences<-cov(differences)

varcov_newton_direction<-cov(newton_direction_values)

#Parts 3 and 4
sample_varcov<-cov(grad_values) #sample variance-covariance of gradient values
sample_avg_hess<-matrix(0,nrow=2,ncol=2) #sample average of Hessian matrices
sample_avg_hess[1,1]<-mean(hessian_values[1,1,])
sample_avg_hess[1,2]<-mean(hessian_values[1,2,])
sample_avg_hess[2,1]<-mean(hessian_values[2,1,])
sample_avg_hess[2,2]<-mean(hessian_values[2,2,])
#Can see that these are quite close
sample_varcov
sample_avg_hess
#Compute inverses
inv_grad_varcov<-solve(sample_varcov)
inv_avg_hess<-solve(sample_avg_hess)
inv_grad_varcov
inv_avg_hess
#Can see that these are very close to each other!
