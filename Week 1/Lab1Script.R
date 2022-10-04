##### Question 1 #####
# Give example of vector theta1/theta2 not working
a<-dat$age[1:17] # get vector of age midpoints
k<-(dat$population[1:17]+dat$population[18:34])/100000 # sizes
theta1=seq(-5,-2,length.out=34) # theta values to test
theta2=seq(0.05,0.2,length.out=34)
nll(theta1,theta2,y,a,k) 

# Need to fix the negative log-likelihood function
nll<-function(theta1,theta2,y,a,k){
  
  lin_pred     <- theta1+theta2*a
  
  mu           <- k*exp(lin_pred)
  
  sum_log_dens <- sum(dpois(x    = y,
                            lambda = mu,
                            log  = TRUE))
  
  -sum_log_dens
}
nll<-Vectorize(FUN=nll,vectorize.args = c("theta1","theta2")) # Vectorize
nll(theta1,theta2,y,a,k) # Make sure that it works

##### Question 2 #####
dat$rate<-dat$deaths_COVID/(dat$population/100000) # Add new variable
par(mfrow=c(1,2))
plot(dat$age[1:17],log(dat$rate[1:17]),xlab='Age',ylab='log(Death rate)',main='Male')
plot(dat$age[18:34],log(dat$rate[18:34]),xlab='Age',ylab='log(Death rate)',main='Female')
fitted<-lm(log(rate)~age,data=dat[1:17,])
fitted$coefficients
