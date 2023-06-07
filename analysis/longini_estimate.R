library(tidyverse)

data <- read_csv("../results/test_results_050323/pptransmission_10inc_results.csv")
data_HH <- read_csv("../results/test_results_050323/pptransmission_10inc_demg.csv")

proportion_HH <- function(sim_data_HH,sim_data,relative=TRUE){
	# household sizes in the data
	k <- sort(unique(sim_data_HH$SIZE))
	# how many households of each size
	n <- table(sim_data_HH$SIZE)
	# matrix with rows=number of infections, cols=household size
	a <- matrix(0,nrow=max(k)+1,ncol=length(k))
	if (length(unique(sim_data$HH))>1){
		cases_per_HH <- table(sim_data$HH,sim_data$SIZE)
	} else {
		cases_per_HH <- table(sim_data$HH)
		cases_per_HH <- t(cases_per_HH)
		colnames(cases_per_HH) <- unique(sim_data$SIZE)[[1]]
	}
	for (idx in 1:length(k)){
		size <- k[idx]
		if (as.character(size) %in% colnames(cases_per_HH)){
			case_counts <- table(cases_per_HH[,as.character(size)])
			for (cases in labels(case_counts)[[1]]){
				a[as.integer(cases)+1,idx] <- case_counts[cases]
			}
		}
	}
	# above loop gives incorrect zero-cases, instead find zero-case
	# households by finding households not included in case data
	a[1,] <- 0
	for (idx in 1:length(k)){
		size <- k[idx]
		for (HH in sim_data_HH[sim_data_HH$SIZE==size,]$HH){
			if (!(HH %in% sim_data$HH)){
				a[1,idx] <- a[1,idx] + 1
			}
		}
	}
	# B_est is probability of escaping infection from community transmission
	B_est <- sum(n*(a[1,]/n)**(1/k))/sum(n)

	# phi is avg number infected per household
	phi <- sum((0:max(k))*rowSums(a))/sum(n)
	# theta is the household attack rate
	theta <- sum((0:max(k))*rowSums(t(t(a)/k)))/sum(n)
	# Q_est is probability of escaping infection from household transmission
	# estimator can give probabilities greater than 1 so set max of 1
	Q_est <- min(1,((1-theta)/B_est)**(1/phi))

	# give relative probability of household infectiohn
	if (relative){
		return((1-Q_est)/(1-B_est))
	}

	# proportion of household infection expressed as conditional probability
	# of household infection given infection at all
	p_HH <- (1-Q_est)/(1-B_est*Q_est)
	return(p_HH)
}

ps <- numeric(length(unique(data$i)))
for (simn in unique(data$i)){
# simn <- 1
	sim_data_HH <- data_HH[data_HH$i==simn,]
	sim_data <- data[data$i==simn,]
	ps[simn] <- proportion_HH(sim_data_HH,sim_data,relative=TRUE)
}
rel_p_hh <- median(ps)
print(rel_p_hh)
