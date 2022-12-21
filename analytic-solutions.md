# Analytical solutions for HEV simulations

## General approach for estimating attack rate
Find the theoretical probability of secondary infection in a susceptible household member in a household with $n$ total individuals given that there is 1 reported infection in a household $i$.  

$$P(\text{secondary infection}) = 1 - P(\text{no secondary infection})$$

$$P(\text{no secondary infection}) = \prod_{t=1}^{38} P(\text{no secondary infection on day }t)$$


### Person-person transmission
Individuals can be infected from the community ( $I_{c,t}$ = number of current infections in the community on day $t$):  

$$risk_{c} = \frac{\beta_{p}I_{c,t}}{(1000-n)}$$

or from household contacts ( $I_{i,t}$ = number of current infections in the household on day $t$ and $r$ = household relative risk):   

$$risk_{h} = \frac{r\beta_{p}I_{i,t}}{n}$$ 

$$P(\text{no secondary infection on day }t) = 1-(risk_{h}+risk_{c}-risk_{h}risk_{c})$$  

$$P(\text{secondary infection}) = 1 - \prod_{t=1}^{38} (1-(risk_{h}+risk_{c}-risk_{h}risk_{c}))$$   
  
  
Ex. Assume 1 infectious person in a household, no prior infections in the household or community - a susceptible household member needs to avoid infection for the approximately 7 day infectious period. Risk of infection is 0 on all other days of the 38 day window.  

$$P(\text{secondary infection}) = 1 - (1-(risk_{h}+risk_{c}-risk_{h}risk_{c}))^{7} = 0.44$$    

### Environmental transmission
Individuals can only be infected through contact with a contaminated environmental source:

$$risk_{e} = \beta_{e}$$  

$$P(\text{no secondary infection}) = \prod_{t=1}^{38} (1-risk_{e}) = (1-risk_{e})^{38}$$ 

$$P(\text{secondary infection}) = 1-(1-risk_{e})^{38}$$

# Ignore!
### Combined person-person and environmental transmission
Individuals can be infected through person-person contact or through environmental transmission:

$$\text{risk of infection} = risk_{h}+risk_{c}+risk_{e}-risk_{h}risk_{c}-risk_{c}risk_{e}-risk_{h}risk_{e}+risk_{h}risk_{c}risk_{e}$$

$$P(\text{no secondary infection in household member }j) = (1-\text{risk of infection})^{7}$$ 
