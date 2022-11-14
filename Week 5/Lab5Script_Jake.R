#Read in the data, 20 obs. of one variable
dat_gamma<-read.table(url("https://people.bath.ac.uk/kai21/ASI/data/gamma_sample.txt"),header = T)

#### Reparametrisation and optimisation ####
#Consider parametrisation given by lambda=(log(alpha),log(beta))
#This takes onto the unconstrained real domain, where we can use optim()

nll_lambda<-function(lambda,data){ #Compute negative log-likelihood
  logdata<-log(data)
  n<-max(dim(data))
  nll<- -n*lambda[2]*exp(lambda[1])+n*lgamma(exp(lambda[1]))-(exp(lambda[1])-1)*sum(logdata)+exp(lambda[2])*sum(data)
  return(nll)
}

nll_lambda_grad<-function(lambda,data){ #Compute gradient of neg. log-likelihood
  logdata<-log(data)
  n<-max(dim(data))
  grad<-matrix(0,nrow=2,ncol=1)
  grad[1]<- -n*lambda[2]*exp(lambda[1])+n*exp(lambda[1])*digamma(exp(lambda[1]))-exp(lambda[1])*sum(logdata)
  grad[2]<- -n*exp(lambda[1])+exp(lambda[2])*sum(data)
  return(grad)
}

#Now use optim to minimise function with respect to unconstrained parameters
fit_lambda<-optim(par     = c(0,0),
                  fn      = nll_lambda,
                  gr      = nll_lambda_grad,
                  data    = dat_gamma,
                  method = "BFGS",
                  hessian = TRUE)
#Gives us the MLE for lambda, which can be transformed back into (alpha,beta)
lambda_star<-fit_lambda$par
shape_rate_star<-exp(lambda_star)

#Can also find Hessian of lambda and using Jacobian, that of (alpha,beta)
#Note that in gamma distribution example, Hessian and fisher information are equal
hess_lambda_observed<-fit_lambda$hessian
J<-diag(exp(-c(fit_lambda$par))) #This is the Jacobian, note the symmetry and remember we are going back from lambda to theta!
hess_theta_observed<-t(J)%*%hess_lambda_observed%*%J 

#Can get inverse of this observed fisher info matrix (in terms of (alpha,beta))
inverse_fisher_theta<-solve(hess_theta_observed) #Has reasonable values

#### Hypothesis test - Confidence intervals ####
#For hypothesis test, use reparametrisation gamma=(alpha-beta^2,beta)
MLE<-c(shape_rate_star[1]-shape_rate_star[2]^2,shape_rate_star[2]) #find reparametrised MLE
J_2<-matrix(0,nrow=2,ncol=2) #Form Jacobian in this case
J_2<-J_2+diag(c(1,1))
J_2[1,2]=2*shape_rate_star[2]
hess_gamma_observed<-t(J_2)%*%hess_theta_observed%*%J_2 #Find reparametrized Fisher info
inverse_fisher_gamma<-solve(hess_gamma_observed)


#We now have all the ingredients required for hypothesis test
lower_bound=MLE[1]-qnorm(0.025,lower.tail = FALSE)*sqrt(inverse_fisher_gamma[1,1])
upper_bound=MLE[1]+qnorm(0.025,lower.tail = FALSE)*sqrt(inverse_fisher_gamma[1,1])
confidence_interval<-c(lower_bound,upper_bound)
confidence_interval
#This contains zero, insufficient evidence to reject null hypothesis

#### Hypothesis test - GLRT ####
#Use invariance property of MLE for unrestrained problem, remember GLRT uses log-likelihood (x-1)
MLE_likelihood_unrestrained<- -nll_lambda(fit_lambda$par,dat_gamma)

#Need to define negative log-likelihood and gradient in restrained case
#Note that we still have to reparametrise to use optim()
#First function is for negative log-likelihood, second is for gradient
nll_restrained<-function(lambda,data){
  logdata<-log(data)
  n<-max(dim(data))
  nll<- -(n/2)*lambda*exp(lambda)+n*lgamma(exp(lambda))-(exp(lambda)-1)*sum(logdata)+exp(lambda/2)*sum(data)
}

nll_restrained_gradient<-function(lambda,data){
  logdata<-log(data)
  n<-max(dim(data))
  grad<- -(n/2)*lambda*exp(lambda)-(n/2)*exp(lambda)+n*exp(lambda)*digamma(exp(lambda))-exp(lambda)*sum(logdata)+0.5*exp(lambda/2)*sum(data)
  return(grad)
}

#Now optimise for restrained MLE and get value of negative log-likelihood at this point
fit_lambda_1<-optim(par   = 0,
                  fn      = nll_restrained,
                  gr      = nll_restrained_gradient,
                  data    = dat_gamma,
                  method = "BFGS")
MLE_likelihood_restrained<- -nll_restrained(fit_lambda_1$par,dat_gamma)

test_statistic<-2*(MLE_likelihood_unrestrained-MLE_likelihood_restrained) #small value of t
critical_value<-qchisq(p=.95,df=1)
test_statistic>critical_value #this is false, we have insufficient evidence to reject null as before