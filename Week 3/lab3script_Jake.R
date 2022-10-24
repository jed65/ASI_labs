library(tidyverse)
covid_districts<- read_csv(file = "https://people.bath.ac.uk/kai21/ASI/data/local_authority_covid_deaths.csv")
covid_districts %>% select(local_authority,death_rate) -> covid_districts #select only the area and number of deaths

#### Question 1 ####
N=10000 #number of samples required
samples<-matrix(0,nrow=N,ncol=311) #matrix to hold samples
for (i in 1:N){ #for loop to generate samples
  samples[i,]<-rnbinom(311,10,0.05) #call rnbinom() to give us samples
} 

#### Question 2 ####
#Need the functions from the notes
#First three functions are for negative log-likelihood
loglikelihood_nbinom2<- function(mean,size,data){
  
  sum_log_dens <- function(x,size,mean){
    sum(dnbinom(x,size,mu=mean,log = TRUE))  
  }
  
  map2_dbl(.x = size,
           .y = mean,
           .f = sum_log_dens,
           x = data)
  
}
loglikelihood_nbinom4<- function(logmean,logsize,data){
  
  mean = exp(logmean)
  size = exp(logsize)
  
  loglikelihood_nbinom2(mean,size,data)
  
}
negloglik_fn<-function(theta = c(0,0), data = 1){
  -loglikelihood_nbinom4(logmean = theta[1], logsize = theta[2],data = data)
}

#Next part is for gradient function (also from the notes)
loglik_expr <- expression(lgamma(y+exp(theta2))-lgamma(exp(theta2))-lgamma(y+1) + exp(theta2)*theta2 + y*theta1 -(exp(theta2)+y)*log(exp(theta1)+exp(theta2)) )
negloglik_deriv <- deriv(expr         = loglik_expr,
                         namevec      = c("theta1","theta2"),
                         function.arg = c("theta1","theta2","y"),
                         hessian      = TRUE)
negloglik_grad <- function(theta = c(0,0),
                           data = 1){
  
  aux  <- negloglik_deriv(theta1 = theta[1],
                          theta2 = theta[2],
                          y      = data)
  
  grad <- apply(attr(aux,"gradient"),2,sum)
  
  -grad
}

#Now we have these, can continue with finding MLEs
parameter_values<-matrix(0,nrow=N,ncol=2) #to store MLEs
grad_values<-matrix(0,nrow=N,ncol=2) #to store gradient at MLEs
hessian_values<-array(0,dim=c(2,2,N)) #store Hessian at MLEs
newton_direction_values<-matrix(0,nrow=N,ncol=2) #store Newton descent directions
for (i in 1:N){
  run<-optim(par=c(0,0),fn=negloglik_fn,gr=negloglik_grad,data=samples[i,],
             method="BFGS",hessian=TRUE) #perform optimisation
  parameter_values[i,]<-run$par #store MLEs
  grad_values[i,]<-negloglik_grad(run$par,samples[i,]) #store gradient at MLE
  hessian_values[,,i]<-run$hessian #store hessian at MLE
  newton_direction_values[i,]<- -solve(run$hessian)%*%grad_values[i,] #calculate newton direction
}
#There is a warning suggesting there could be some NaNs but 
#I checked using any(is.nan()) and got false for all 4 storage arrays


