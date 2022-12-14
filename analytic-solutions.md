# Analytical solutions for HEV simulations

## General approach
Find the theoretical probability of any secondary infection among susceptible household members $j=1...$ in a household with n total members given that there is 1 reported infection in household $i$ at time $t$.  

$P(any secondary infection) = 1 - P(no secondary infection)$
$P(no secondary infection) = \prod_{j = 1}^{S_{i,t}} P(no secondary infection in household member i)$
$P(no secondary infection in household member j) = P(household member j not infected by any possible source such that infection would begin anywhere between t+7 to t+45 days)$

If someone were currently infectious with person-person transmission - they can be infected from the community ($risk_{c} = household rr*beta_{pp}*I_{c,t}/(1000-n)$) or from household contacts ($risk_{h} = household rr*beta_{pp}*I_{i,t}/(n-1)$)  
$P(no secondary infection) = (1-(risk_{h}+risk_{c}-risk_{h}risk_{c}))**7