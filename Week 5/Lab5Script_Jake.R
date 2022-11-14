#Read in the data, 20 obs. of one variable
dat_gamma<-read.table(url("https://people.bath.ac.uk/kai21/ASI/data/gamma_sample.txt"),header = T)

#Explicit functions for likelihood and gradient with theta=(alpha,beta) (shape and rate parameters)
nll_fn<-function(theta=c(0,0),data=1){
  logdata<-log(data)
  nll<- -20*theta[1]*log(theta[2])+20*lgamma(theta[1])-(theta[1]-1)*sum(logdata)+theta[2]*sum(data)
  return(nll)
}

nll_grad<-function(theta=c(0,0),data=1){
  logdata<-log(data)
  nll_grad_theta1<- -20*log(theta[2])+20*digamma(theta[1])-sum(logdata)
  nll_grad_theta2<- -20*(theta[1]/theta[2])+sum(data)
  grad_vector<-matrix(0,nrow=2,ncol=1)
  grad_vector[1]=nll_grad_theta1
  grad_vector[2]=nll_grad_theta2
  return(grad_vector)
}

#First we calculate an asymptotic confidence interval
#Now use optim to fit to the data
fit_alpha_beta<-optim(par=c(0.5,0.5),fn=nll_fn,gr=nll_grad,
                      data=dat_gamma,method="BFGS",hessian=TRUE)

fisher_shape_rate_observed<-fit_alpha_beta$hessian #use observed information
inverse_fisher_shape_rate_observed<-solve(fisher_shape_rate_observed) #find inverse

#Now consider reparametrisation given by g_(1)^(-1)(alpha,beta)=alpha-beta^2 and g_(1)^(-1)(alpha,beta)=beta
#MLE for alpha-beta^2
lambda_1_star<-fit_alpha_beta$par[1]-fit_alpha_beta$par[2]^2

#Inverse fisher info matrix
J_g_inverse=matrix(0,nrow=2,ncol=2)
J_g_inverse[1,1]=1; J_g_inverse[2,2]=1; J_g_inverse[1,2]=-2*fit_alpha_beta$par[2]
#Find reparametrised inverse fisher matrix
inverse_fisher_reparametrised<-J_g_inverse%*%inverse_fisher_shape_rate_observed%*%t(J_g_inverse)

#Finally, find the confidence interval
lower_bound<-lambda_1_star-qnorm(0.025,lower.tail = FALSE)*sqrt(inverse_fisher_reparametrised[1,1])
upper_bound<-lambda_1_star+qnorm(0.025,lower.tail = FALSE)*sqrt(inverse_fisher_reparametrised[1,1])


#Print to screen
lower_bound
upper_bound
#Gives the interval (-14.60623,8.120496), which contains 0 so evidence supporting null hypothesis


#Now we need to consider the generalised likelihood ratio test