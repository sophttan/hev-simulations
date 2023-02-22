# British Boarding School SIR, v1
# Date- 7/18/16


# Load R packages  
require(bbmle) # General MLE
require(deSolve) # Differential equation solver 

# Test data v3- British boarding school flu case data (N=764)
incidence <- c(3,6,25,73,222,294,258,237,191,125,69,27,11,4)
week <- seq(1,14)
data<- data.frame(week, incidence)
plot(data$week, data$incidence, type='b',col='red')

# Simple SIR model 
sir.model.closed <- function (t, y, params) { # Note: this input format is required for ode function (used in "prediction" function)
  # x- initial guesses 
  S <- y[1]
  I <- y[2]
  R <- y[3]
  with(
    as.list(params),
    {
      dS <- -beta*S*I
      dI <- beta*S*I-gamma*I
      dR <- gamma*I
      dx<- c(dS, dI, dR)
      list(dx)
    }
  )
}

prediction <- function (params, times) {
  xstart <- params[c("S.0","I.0", "R.0")]
  # ODE function
  # y- initial conditions
  # func- set of ODEs
  # params- passed to function 
  out <- ode(
    func=sir.model.closed,
    y=xstart,
    times= c(0, times),
    parms=params,
    hmax=1/120
  )
  out[-1,3] # return the I variable only
}

# Define likelihood function with poisson distribution (count data)
poisson.loglik <- function (params, data) {
  times <- data$week # convert to time-scale year
  pred <- prediction(params,times)
  #sum(dpois(x=data$incidence,lambda=pred[-1],log=TRUE))
  poiss <- sum(dpois(x=data$incidence,lambda=pred*params["N"],log=TRUE)) # include "p" for reporting probability 
  print(poiss)
  poiss
  }

params <- c(S.0=763/764,I.0=1/764, R.0=0/764, gamma=NA,beta=NA, N=764)
 #logit <- function (p) log(p/(1-p))      # the logit transform
 #ilogit <- function (x) 1/(1+exp(-x))    # inverse logit

# f
# Objective: How does a set of model parameters perform for poisson negative log-likelihood?
likelihood <- function (par_initial) {
  par <- params
  par[c("gamma")] <- par_initial[1]
  par[c("beta")] <- par_initial[2]
  #print(par)
  poisson.loglik(par,data)
}







# Initial guess
# Make sure to vary this initial guess to determine if guess robust. 
#guess <- list(log.beta=log(0.005),logit.p=logit(0.2))
#guess <- list(gamma=0.45,beta=1.5)

# fit0 <- mle2(f,start=guess)

# Metroplalis 
library(mcmc)

#mcmc <- metrop(f, initial=c(0.45,1.5), nbatch=100)

mcmc <- metrop(likelihood, initial=c(0.1,1), nbatch=100)

params["beta"] <- mcmc$final[2]
params["gamma"] <- mcmc$final[1]
times <- c(data$week)
model.pred <- prediction(params,times)*params["N"]
plot(incidence~week,data=data,type='b',col='red') 
lines(times,model.pred,type='l')


# Metropalis-Hasting
library(MHadaptive)


mcmc_MH <- Metro_Hastings(li_func=likelihood, pars=c(0.1,1), iterations=1000, burn_in=100)



# MCMCpack 
post_sample <- MCMCmetrop1R(fun = function(par) { likelihood(par) }, theta.init= c(0.1,1), mcmc = 500, burnin = 100)








## MCMC- MC self code


# Prior distribution
prior <- function(par_initial){
  gamma = par_initial[1]
  beta = par_initial[2]

  gamma_prior = dnorm(gamma, mean=0.5, sd=0.1, log = T)
  beta_prior = dnorm(beta, mean=1.5, sd = 0.25, log = T)
  return(gamma_prior+beta_prior)
}


posterior <- function(param){
  return (likelihood(param) + prior(param))
}


proposalfunction <- function(param){
  return(rnorm(2,mean = param, sd= c(0.1,0.1)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

sim_size <- 3000
startvalue = c(0.1, 1)
chain = run_metropolis_MCMC(startvalue, sim_size)

burnIn = 1000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))



params["beta"] <- chain[sim_size, 2]
params["gamma"] <- chain[sim_size, 1]
times <- c(data$week)
model.pred <- prediction(params,times)*params["N"]
plot(incidence~week,data=data,type='b',col='red') 
lines(times,model.pred,type='l')

