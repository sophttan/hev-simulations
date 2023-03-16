library(tidyverse)

time = 365
inc = 28
inf = 7
pop = 1000 # population

create_hh <- function() {
  # randomly sample household sizes such that total population is 1000 individuals
  hh_size <- sample(x = c(3, 4, 5, 6), size=340, replace=T)
  # keep households such that total population is < 1000
  hh_size <- hh_size[which(cumsum(hh_size) < pop)]
  
  leftover <- pop-sum(hh_size)
  if(leftover < 3) {
    hh <- 1:length(hh_size)
    sampled <- sample(hh[hh_size<6], leftover)
    hh_size[sampled] <- hh_size[sampled] + 1
  } else {
    hh_size <- c(hh_size, leftover)
  }
  
  return(hh_size)
}

SEIR <- function(beta_H, beta_C, inc, inf, verbose = 0) {
  hh_size = create_hh()
  
  # ID: ID of individual
  # HHsize: size of individual's household
  # HH: ID of individual's household
  # S: susceptibility status
  # E: exposed status
  # Ecounter: number of days since exposed
  # Icounter: infectious status
  # I_count: number of days since infectious
  # R: recovered status
  # inc: incubation period
  # inf: infectious period
  
  d = data.frame(
    ID=1:pop, 
    HHsize=rep(hh_size, times=hh_size),
    HH=rep(1:length(hh_size), times=hh_size), 
    S=c(0, rep(1, pop-1)), 
    E=c(1, rep(0, pop-1)),
    Ecounter=c(1, rep(0, pop-1)), 
    I=0,
    I_count=0, 
    R=0, 
    inc=c(round(rnorm(1, inc, 2)), rep(0, pop-1)),
    inf=0
  )
  
  # create frame for storing results
  results <- d[,1:3] %>% mutate(type=NA, time=NA)
  results$type[1] = 0
  results$time[1] = 0
  
  for(t in 1:time) { # (?) what does time here represent?
    if(verbose) {
      if(t %% 10 == 0) {
        paste0(t) # (?) progress bar?
      }
    }
    
    # Anyone who has been infectious for as many days as their infectious
    # period is now recovered.
    recovered <- d$inf>0 & d$Icounter==d$inf
    if(sum(recovered,na.rm=T)>0) {
      d$R[recovered] <- 1
      d$I[recovered] <- 0
      d$Icounter[recovered] <- 0 
    }
    
    # Anyone who has been incubating for as many days as their incubation
    # period is now infectious.
    new_inf <- d$inc>0 &d$Ecounter==d$inc
    if(sum(new_inf,na.rm=T)>0) {
      random_inf <- rnorm(sum(new_inf, na.rm=T), mean=inf, sd=1) %>% round()
      d$I[new_inf] <- 1
      d$inf[new_inf] <- random_inf
      d$E[new_inf] <- 0
      d$Ecounter[new_inf] <- 0 
    }
    
    # I_H is the number of infections in each household.
    # I_C is the number of infections outside a given household.
    I_H = d %>% group_by(HH) %>% summarise(sum = sum(I, na.rm = T))
    summ = data.frame(I_H=I_H, 
                      I_C = sum(d$I) - I_H
    )
    
    # dd is a frame where each individual is assigned their household's
    # I_H and I_C numbers.
    dd = d %>% select(HH, S)
    dd = dd %>% mutate(I_H = summ$I_H.sum[HH])
    dd = dd %>% mutate(I_C = summ$I_C.sum[HH])
    
    # Calculate household risk and community risk
    risk_H <- beta_H*d$S*dd$I_H / pop
    risk_C <- beta_C*d$S*dd$I_H / pop
    
    new_inf_H <- rbinom(nrow(d), 1, risk_H)
    new_inf_C <- rbinom(nrow(d), 1, risk_C)
    new_exposed = (new_inf_H == 1) | (new_inf_C == 1)
    
    if(sum(new_exposed)>0) {
      d$E[new_exposed==1] <- 1
      d$inc[new_exposed==1] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
      results <- results %>% mutate(time = ifelse(new_exposed==1, t, time))
      
      results$type[!((new_inf_H == 1) & (new_inf_C == 1)) | !is.na(results$type)] <- 'B'
      results$type[!(new_inf_H == 1) | !is.na(results$type)] <- 'H'
      results$type[!(new_inf_H == 1) | !is.na(results$type)] <- 'C'
      results$time[!(new_exposed == 1) | !is.na(results$time)] <- t
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
  }
  return(results)
}

# for testing SEIR
# beta <- 0.0011
# sims <- 1000
# num_hh <- rep(0, sims)
# inc <- rep(0, sims)
# inf <- rep(0, sims)
# 
# SEIR(beta, beta, inc, inf)

metrics = function(results) {
  state = results %>% slice(which(!is.na(results$time)))
  idc = nrow(state)/1000
  
  sar = NA
  f = function(x) {
    !all(is.na(x))
  }
  if (idc != 0) {
    num_primary = sum(results %>% group_by(HH) %>% summarise(time_sum = sum(time, na.rm = T)) %>% ungroup() %>% select(time_sum) > 0) # households that were infected
    idx = results %>% select(type) %>% sapply(f)
    num_contact = results %>%
      group_by(HH) %>%
      summarise(size = sum(HHsize, na.rm = T)) %>%
      ungroup() %>%
      select(size) %>%
      slice(which(idx)) %>%
      mutate(size2 = size^(1/2)) %>%
      sum() # total people in all those households
    sar = dim(state[(state$type == 'H') | (state$type == 'B')])[1] / (num_contact - num_primary)
  }
  return(c(idc, sar))
}

# for the next function to work we need the metric function to output
# something that can be converted into a matrix

score = function(results, target) {
  return(sum(metrics(results) - target)^2)
}

# Create likelihood from the score of the state
likelihood = function(state) {
  beta_H = state[1]
  beta_C = state[2]
  
  results = SEIR(beta_H, beta_C, inc, inf)
  liks = -log(score(results, target))
  # (?) different from the original python code which seems to have an unnecessary for loop
  return(liks)
}

#### Metropolis algorithm ####

# Proposal function
q = function(state) {
  beta_H = state[1]
  beta_C = state[2]
  r = c(beta_H^2, beta_C^2 / 0.0001)
  v = c(beta_H^2, beta_C^2 / 0.0001)
  
  return(pmax(rgamma(n = 2, shape = r, scale = 1 / v), 1e-3))
  # (?) the parameters of the above gamma distribution may be wrong, but i
  # assumed r to be the shape parameter
}

# MCMC
metropolis = function(start, num_iter) {
  chain = matrix(0, num_iter + 1, 2)
  liks = matrix(0, num_iter + 1, 2)
  
  # Initialize current state.
  curr = start
  curr_lik = likelihood(curr)
  
  # Initialize best state.
  best = curr
  best_lik = curr_lik
  for(i in 1:num_iter) {
    # Save the current state and its likelihood.
    chain[i,] = curr
    liks[i,] = curr_lik
    
    # Print current state and likelihood.
    paste0(i, '\t', curr, " ", signif(curr_lik,3), '\t')
    
    # Get a proposed state and calculate its likelihood.
    prop = q(curr)
    prop_lik = likelihood(prop)
    paste0(prop, " ", signif(prop_lik,3), '\t')
    
    # Compute the ratio of the scores of the two states
    # and flip a coin.
    r = exp(prop_lik - curr_lik)
    p = runif(1)
    paste(signif(r, 3), signif(p, 3))
    library(plyr)
    library(tidyverse)
    
    time = 365
    inc = 28
    inf = 7
    pop = 1000 # population
    
    create_hh <- function() {
      # randomly sample household sizes such that total population is 1000 individuals
      hh_size <- sample(x = c(3, 4, 5, 6), size=340, replace=T)
      # keep households such that total population is < 1000
      hh_size <- hh_size[which(cumsum(hh_size) < pop)]
      
      leftover <- pop-sum(hh_size)
      if(leftover < 3) {
        hh <- 1:length(hh_size)
        sampled <- sample(hh[hh_size<6], leftover)
        hh_size[sampled] <- hh_size[sampled] + 1
      } else {
        hh_size <- c(hh_size, leftover)
      }
      
      return(hh_size)
    }
    
    SEIR <- function(beta_H, beta_C, inc, inf, verbose = 0) {
      hh_size = create_hh()
      
      # ID: ID of individual
      # HHsize: size of individual's household
      # HH: ID of individual's household
      # S: susceptibility status
      # E: exposed status
      # Ecounter: number of days since exposed
      # Icounter: infectious status
      # I_count: number of days since infectious
      # R: recovered status
      # inc: incubation period
      # inf: infectious period
      
      d = data.frame(
        ID=1:pop, 
        HHsize=rep(hh_size, times=hh_size),
        HH=rep(1:length(hh_size), times=hh_size), 
        S=c(0, rep(1, pop-1)), 
        E=c(1, rep(0, pop-1)),
        Ecounter=c(1, rep(0, pop-1)), 
        I=0,
        I_count=0, 
        R=0, 
        inc=c(round(rnorm(1, inc, 2)), rep(0, pop-1)),
        inf=0
      )
      
      # create frame for storing results
      results <- d[,1:3] %>% mutate(type=NA, time=NA)
      results$type[1] = 0
      results$time[1] = 0
      
      for(t in 1:time) { # (?) what does time here represent?
        if(verbose) {
          if(t %% 10 == 0) {
            paste0(t) # (?) progress bar?
          }
        }
        
        # Anyone who has been infectious for as many days as their infectious
        # period is now recovered.
        recovered <- d$inf>0 & d$Icounter==d$inf
        if(sum(recovered,na.rm=T)>0) {
          d$R[recovered] <- 1
          d$I[recovered] <- 0
          d$Icounter[recovered] <- 0 
        }
        
        # Anyone who has been incubating for as many days as their incubation
        # period is now infectious.
        new_inf <- d$inc>0 &d$Ecounter==d$inc
        if(sum(new_inf,na.rm=T)>0) {
          random_inf <- rnorm(sum(new_inf, na.rm=T), mean=inf, sd=1) %>% round()
          d$I[new_inf] <- 1
          d$inf[new_inf] <- random_inf
          d$E[new_inf] <- 0
          d$Ecounter[new_inf] <- 0 
        }
        
        # I_H is the number of infections in each household.
        # I_C is the number of infections outside a given household.
        I_H = d %>% group_by(HH) %>% summarise(sum = sum(I, na.rm = T))
        summ = data.frame(I_H=I_H, 
                          I_C = sum(d$I) - I_H
        )
        
        # dd is a frame where each individual is assigned their household's
        # I_H and I_C numbers.
        dd = d %>% select(HH, S)
        dd = dd %>% mutate(I_H = summ$I_H.sum[HH])
        dd = dd %>% mutate(I_C = summ$I_C.sum[HH])
        
        # Calculate household risk and community risk
        risk_H <- beta_H*d$S*dd$I_H / pop
        risk_C <- beta_C*d$S*dd$I_H / pop
        
        new_inf_H <- rbinom(nrow(d), 1, risk_H)
        new_inf_C <- rbinom(nrow(d), 1, risk_C)
        new_exposed = (new_inf_H == 1) | (new_inf_C == 1)
        
        if(sum(new_exposed)>0) {
          d$E[new_exposed==1] <- 1
          d$inc[new_exposed==1] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
          results <- results %>% mutate(time = ifelse(new_exposed==1, t, time))
          
          results$type[!((new_inf_H == 1) & (new_inf_C == 1)) | !is.na(results$type)] <- 'B'
          results$type[!(new_inf_H == 1) | !is.na(results$type)] <- 'H'
          results$type[!(new_inf_H == 1) | !is.na(results$type)] <- 'C'
          results$time[!(new_exposed == 1) | !is.na(results$time)] <- t
        }
        
        d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
        d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
        d$S[d$E==1] <- 0
      }
      return(results)
    }
    
    # for testing SEIR
    # beta <- 0.0011
    # sims <- 1000
    # num_hh <- rep(0, sims)
    # inc <- rep(0, sims)
    # inf <- rep(0, sims)
    # 
    # SEIR(beta, beta, inc, inf)
    
    metrics = function(results) {
      state = results %>% slice(which(!is.na(results$time)))
      idc = nrow(state)/1000
      
      sar = NA
      f = function(x) {
        !all(is.na(x))
      }
      if (idc != 0) {
        num_primary = sum(results %>% group_by(HH) %>% summarise(time_sum = sum(time, na.rm = T)) %>% ungroup() %>% select(time_sum) > 0) # households that were infected
        idx = results %>% select(type) %>% sapply(f)
        num_contact = results %>%
          group_by(HH) %>%
          summarise(size = sum(HHsize, na.rm = T)) %>%
          ungroup() %>%
          select(size) %>%
          slice(which(idx)) %>%
          mutate(size2 = size^(1/2)) %>%
          sum() # total people in all those households
        sar = dim(state[(state$type == 'H') | (state$type == 'B')])[1] / (num_contact - num_primary)
      }
      return(c(idc, sar))
    }
    
    # for the next function to work we need the metric function to output
    # something that can be converted into a matrix
    
    score = function(results, target) {
      return(sum(metrics(results) - target)^2)
    }
    
    # Create likelihood from the score of the state
    likelihood = function(state) {
      beta_H = state[1]
      beta_C = state[2]
      #avg(metrics(results))
      #then estimate score relative to the target
      #then estimate likelihood
      
      Nresults = SEIR(beta_H, beta_C, inc, inf)
      Nresults = Nresults %>% mutate(df = 1)
      
      for(i in 2:N) {
        results = SEIR(beta_H, beta_C, inc, inf)
        results = results %>% mutate(df = i)
        Nresults = cbind(Nresults, results)
      }
      
      split_results <- split(Nresults, f = Nresults$df)
      avg_results <- Reduce(`+`, mget(paste0("split_results$", 1:N)))/N
      liks[i] = -log(score(avg_results, target))
      return(liks)
    }
    
    #### Metropolis algorithm ####
    
    # Proposal function
    q = function(state) {
      beta_H = state[1]
      beta_C = state[2]
      r = c(beta_H^2, beta_C^2 / 0.0001)
      v = c(beta_H^2, beta_C^2 / 0.0001)
      
      return(pmax(rgamma(n = 2, shape = r, scale = 1 / v), 1e-3))
      # (?) the parameters of the above gamma distribution may be wrong, but i
      # assumed r to be the shape parameter
    }
    
    # MCMC
    metropolis = function(start, num_iter) {
      chain = matrix(0, num_iter + 1, 2)
      liks = matrix(0, num_iter + 1, 2)
      
      # Initialize current state.
      curr = start
      curr_lik = likelihood(curr)
      
      # Initialize best state.
      best = curr
      best_lik = curr_lik
      for(i in 1:num_iter) {
        # Save the current state and its likelihood.
        chain[i,] = curr
        liks[i,] = curr_lik
        
        # Print current state and likelihood.
        paste0(i, '\t', curr, " ", signif(curr_lik,3), '\t')
        
        # Get a proposed state and calculate its likelihood.
        prop = q(curr)
        prop_lik = likelihood(prop)
        paste0(prop, " ", signif(prop_lik,3), '\t')
        
        # Compute the ratio of the scores of the two states
        # and flip a coin.
        r = exp(prop_lik - curr_lik)
        p = runif(1)
        paste(signif(r, 3), signif(p, 3))
        
        # Transition if the proposed state is better or
        # if the coin flip succeeds.
        if(p < r) { # (?) use a bernoulli maybe?
          curr = prop
          curr_lik = prop_lik
          
          # If the new likelihood is better than the
          # best we've seen so far, replace the best.
          if(curr_lik > best_lik) {
            best = curr
            best_lik = curr_lik
          }
        }
        
        # Save the chain, best state, and likelihoods
        # so far.
        save(chain, file = 'chain.npy')
        save(liks, file = 'liks.npy')
        save(best, file = 'best.npy')
      }
      return(list(chain, liks, best))
    }
    
    # Solve for optimal values via MCMC
    target = c(0.3, 0.25) # target values
    N = 10 # number of times over which to average likelihood
    
    metropolis_results = metropolis(c(30, 0.12), 10)
    chain = metropolis_results[[1]]
    liks = metropolis_results[[2]]
    best = metropolis_results[[3]]
    save(chain, file = "chain.Rdata")
    save(liks, file = "liks.Rdata")
    save(best, file = "best.Rds")
    
    # Transition if the proposed state is better or
    # if the coin flip succeeds.
    if(p < r) { # (?) wouldn't
      curr = prop
      curr_lik = prop_lik
      
      # If the new likelihood is better than the
      # best we've seen so far, replace the best.
      if(curr_lik > best_lik) {
        best = curr
        best_lik = curr_lik
      }
    }
    
    # Save the chain, best state, and likelihoods
    # so far.
    save(chain, file = 'chain.npy')
    save(liks, file = 'liks.npy')
    save(best, file = 'best.npy')
  }
  return(list(chain, liks, best))
}

# Solve for optimal values via MCMC
target = c(0.3, 0.25) # target values
N = 10 # number of times over which to average likelihood

metropolis_results = metropolis(c(30, 0.12), 10)
chain = metropolis_results[[1]]
liks = metropolis_results[[2]]
best = metropolis_results[[3]]
save(chain, file = "chain.Rdata")
save(liks, file = "liks.Rdata")
save(best, file = "best.Rds")
