source('Schisto_Simulations_Functions.R')

args <- commandArgs(trailingOnly = TRUE)
baseline_prevalence <- as.numeric(args[1])/100
con <- file(paste("seeds",args[1],"_upcap.txt",sep=''))
good_seeds_wNA <- scan(con,what=numeric(),sep=",")
close(con)
good_seeds <- na.omit(good_seeds_wNA)
drug_strategy <- args[2]
coverage <- as.numeric(args[3])/100
static_force <- as.numeric(args[4])/100
binomial_treatment <- as.logical(args[5])
cat(binomial_treatment)
n <- as.numeric(args[6])

seed <- good_seeds[n]
set.seed(seed)

if (baseline_prevalence==0.15){
  lambda_all <- 0.0013
  dispersion_param <- 2.45
  mu_param <- -0.44
} else if (baseline_prevalence==0.3){
  lambda_all <- 0.0025
  dispersion_param <- 2.55
  mu_param <- 0.26
} else if (baseline_prevalence==0.5){
  lambda_all <- 0.0062
  dispersion_param <- 2.48
  mu_param <- 0.8
}

if (binomial_treatment){
  eff1 <- 0.945
  eff2 <- 0.21
} else {
  eff1 <- 0.935
  eff2 <- 0.38
}

if (drug_strategy=='single'){
  treat_weeks <- c(1)
  eff_a <- eff1
  eff_j <- 0
} else if (drug_strategy=='novel'){
  treat_weeks <- c(1)
  eff_a <- eff1
  eff_j <- eff1
} else if (drug_strategy=='double'){
  treat_weeks <- c(1,7)
  eff_a <- eff1
  eff_j <- 0
}
eff_d <- eff2

prevalence_state <- readRDS(paste("baselines/prev",as.character(baseline_prevalence*100),"_upcap/",as.character(seed),".RData",sep=''))
analysis <- strategy_sim(prevalence_state,lambda_all,treat_weeks,eff_a,eff_j,11,delay_efficacy=eff_d,N=500,coverage=coverage,static_force=static_force,binomial_treatment=binomial_treatment)
outcomes <- outcomes_row(analysis,drug_strategy,coverage,prevalence_state)

cat(paste(outcomes,collapse=","),'\n',sep='')
