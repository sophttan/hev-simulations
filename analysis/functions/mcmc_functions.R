# Functions for SEIR model calibration with the Metropolis Algorithm
# Authors: Sophia Tan, Shahzar Rizvi, Nila Cebu

metrics <- function(results) {
  # Incidence is the proportion of the population that became infected.
  idc <- mean(!is.na(results$TIME))
  
  # If incidence is 0, the SAR is undefined.
  sar <- NA
  if (idc != 0) {
    # The SAR is the average SAR for each individual that was infectious.
    sar <- mean(results$I_num / results$S_num, na.rm = T)
  }
  return(c(idc, sar))
}

score <- function(fit, target) {
  # The score is the L² distance of the observed values from the target.
  return(sum((fit - target)^2))
}

# The likelihood is calculated by first averaging the incidence and SAR over n
# simulations with the state parameters. The likelihood is the negative log
# score of the average incidence and SAR.
inf <- 7
likelihood <- function(state, target, n = 300) {
  # If either parameter is nonpositive, do not transition to that state.
  if (any(state <= 0)) {
    return(-Inf)
  }
  # Otherwise, find the average incidence and SAR and compute likelihood.
  vals <- foreach (i = 1:n, .combine = c) %dopar% {
    results <- SEIR(state, inf)
    metrics(results)
  }
  vals <- matrix(vals, n, byrow = T)
  fit <- colMeans(vals)
  return(-log(score(fit, target)))
}

# Proposal function
q <- function(state, sds = c(0.5, 0.005)) {
  # Sample from a multivariate normal distributions centered at the current 
  # state. The SDs roughly correspond to the step-size of the chain for each 
  # parameter.
  return(rnorm(n = 2, mean = state, sd = sds))
}

# MCMC
metropolis <- function(start, target, num_sim, num_iter) {
  path <- matrix(NA, num_iter + 1, 2)
  liks <- rep(NA, num_iter + 1)
  
  # Initialize current state.
  curr <- start
  curr_lik <- likelihood(curr, target, num_sim)
  
  # Initialize best state.
  best <- curr
  best_lik <- curr_lik
  for (i in 1:num_iter) {
    # Save the current state and its likelihood.
    path[i, ] <- curr
    liks[i] <- curr_lik
    
    # Print the current state and its likelihood.
    cat('[', 
        format(curr[1], digits = 5, nsmall = 3), ' ', 
        format(curr[2], digits = 3, nsmall = 4), ']\t', 
        format(curr_lik, digits = 3, nsmall = 3), '\t', sep = '')
      
    # Get a proposed state and calculate its likelihood.
    prop <- q(curr)
    prop_lik <- likelihood(prop, target, num_sim)
    
    # Print the proposed state and its likelihood.
    cat(i, '\t[', 
        format(prop[1], digits = 5, nsmall = 3), ' ', 
        format(prop[2], digits = 3, nsmall = 4), ']\t', 
        format(prop_lik, digits = 3, nsmall = 3), '\t', sep = '')
      
    # Compute the ratio of the scores of the two states and generate a uniform 
    # bit.
    r <- exp(prop_lik - curr_lik)
    p <- runif(1)
    
    # Print the rejection process variables.
    cat(format(r, digits = 3, nsmall = 3), '\t', 
        format(p, digits = 3, nsmall = 3), '\n', sep = '')
    
    # Transition if the proposed state is better or if the coin flip succeeds.
    if (p < r) { 
      curr <- prop
      curr_lik <- prop_lik
      
      # If the new likelihood is better than the best we've seen so far, replace 
      # the best.
      if (curr_lik > best_lik) {
        best <- curr
        best_lik <- curr_lik
      }
    }
    
    # # Save the path, best state, and likelihoods so far.
    # write.table(path, file = paste0('path',  '.txt'), row.names = F, col.names = F)
    # write.table(liks, file = 'liks.txt', row.names = F, col.names = F)
    # write.table(best, file = 'best.txt', row.names = F, col.names = F)
  }
  path[num_iter + 1, ] <- curr
  liks[num_iter + 1] <- curr_lik
  return(list(path, liks, best))
}
