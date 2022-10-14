nll <- function(theta,y,a,k,delta){ # Function to compute neg. log-likelihood
  
  lin_pred     <- theta[1]+theta[2]*a+theta[3]*delta
  
  mu           <- k*exp(lin_pred)
  
  sum_log_dens <- sum(dpois(x = y,
                            lambda = mu,
                            log  = TRUE))
  -sum_log_dens
}

nll_grad <- function(theta,y,a,k,delta){ # Function to compute gradient
  lin_pred     <- theta[1]+theta[2]*a+theta[3]*delta
  mu           <- k*exp(lin_pred)
  
  grad_1 <- sum(mu-y)
  grad_2 <- sum(a*(mu-y))
  grad_3<- sum(delta*(mu-y))
  c(grad_1,grad_2,grad_3)
}

deaths<-read.csv("https://people.bath.ac.uk/kai21/ASI/data/COVID19_MARCH_JUNE.csv",header =TRUE)

deaths$rate<-deaths$deaths_COVID/(deaths$population/100000) #make variable for rate
deaths$sex_logic<-c(rep(1,17),rep(0,17)) #converting sex to a logical operator
deaths$rate_log<-log(deaths$rate) #taking the log of the death rate

multi.fit = lm(rate_log~age+sex_logic, data=deaths) #fit linear model
initial_guess = multi.fit$coefficients #extract starting theta
names(initial_guess) <- NULL #remove names, gives vector theta_0
cat('The starting value from fitting linear model is:', initial_guess, '\n')

# Now to implement the BFGS algorithm
BFGS<-function(theta_0,func,grad,B,
               control=list(maxit=1000,abstol=1e-6,c1=1e-4,c2=0.9),
               y,a,k_in,delta){
  norm_vec<- function(x) sqrt(sum(x^2)) # Compute l2 norm
  k<-1 # Iteration counter
  convergence_indicator <- 0 # 0 if successful, 1 otherwise
  d <- length(theta_0) # Dimension of our problem
  theta_current<-theta_0 # Set current value of theta
  nll_current<-func(theta_current,y,a,k_in,delta) # Evaluate l()
  grad_nll_current<-grad(theta_current,y,a,k_in,delta) # Evaluate grad l()
  theta_iterates<-matrix(0,nrow=control$maxit,ncol=d) # Matrix to hold iterates
  theta_iterates[k,]<-theta_current # Set first line in matrix to initial guess
  k_backtracking <- rep(0,control$maxit) # Keep track of how many backtracks at each iteration
  
  # Now we proceed with the main part of the algorithm
  while (norm_vec(grad_nll_current)>control$abstol){ # Stop criterion
    
    Delta <- -B %*% grad_nll_current # Computes descent direction
    
    # Backtracking section
    alpha<-1 # This is alpha_k in notes
    alpha_min<-0 
    alpha_max<-Inf
    
    # Find values for this value of step length
    theta_proposed<-theta_current+alpha*Delta
    nll_proposed<-func(theta_proposed,y,a,k_in,delta)
    grad_nll_proposed<-grad(theta_proposed,y,a,k_in,delta)
    
    k_bt<-0 # Keeps track of number of backtracks
    # Write out the conditions for modified backtracking
    cond1 <- nll_proposed >= nll_current + control$c1*alpha*crossprod(grad_nll_current,Delta)
    cond2 <- crossprod(grad_nll_proposed,Delta) < control$c2*crossprod(grad_nll_current,Delta)
    
    while(cond1|cond2){ # While either of these conditions holds
      if(cond1){
        alpha_max=alpha
        alpha=(alpha_min+alpha_max)/2
      }else if(cond2){
        alpha_min=alpha
        if(alpha_max==Inf){
          alpha=2*alpha
        }else{
          alpha=(alpha_min+alpha_max)/2
        }
      } 
      # Re-evaluate proposed values for new alpha
      theta_proposed<-theta_current+alpha*Delta
      nll_proposed<-func(theta_proposed,y,a,k_in,delta)
      grad_nll_proposed<-grad(theta_proposed,y,a,k_in,delta)
      # Re-evaluate conditions to go back into loop
      cond1 <- nll_proposed >= nll_current + control$c1*alpha*crossprod(grad_nll_current,Delta)
      cond2 <- crossprod(grad_nll_proposed,Delta) < control$c2*crossprod(grad_nll_current,Delta)
      
      k_bt<-k_bt+1 # Add 1 to number of backtracks
    }
    
    k_backtracking[k]<-k_bt # Store number of backtracks
    # Update and store values
    grad_nll_old<-grad_nll_current
    theta_current<-theta_proposed
    nll_current<-nll_proposed
    grad_nll_current<-grad(theta_current,y,a,k_in,delta)
    # Things we need to update B
    eta<-grad_nll_current-grad_nll_old
    rho <- 1/as.numeric(crossprod(eta,Delta))
    # Update B
    B=(diag(d)-rho*tcrossprod(Delta,eta))%*%B%*%(diag(d)-rho*tcrossprod(eta,Delta))+alpha*rho*tcrossprod(Delta,Delta)
    
    k <- k+1 # Add one to number of iterations
    theta_iterates[k,]<-theta_current # Store the iteration
    
    if (k==control$maxit){
      convergence_indicator<-1
    }
  }
  # Pass back the outputs
  return(list(theta=theta_current,
              nll_value=nll_current,
              convergence=convergence_indicator,
              no_iterations=k,
              iterates=theta_iterates[1:k,],
              inv_hess=B,
              no_backtracks=k_backtracking[1:k]))
}

# Set input arguments from dataset
y = deaths$deaths_COVID
a = deaths$age
delta = deaths$sex_logic
k = deaths$population/100000

# Need to pass an initial matrix B, set this as identity matrix
# There is a default in control so no need to pass
run_BFGS<-BFGS(theta_0=initial_guess,func=nll,grad=nll_grad,B=diag(3),y=y,a=a,k_in=k,delta=delta)
run_BFGS # Display all the information
theta_final<-as.numeric(run_BFGS$theta) # This is the final value of theta from the algorithm


