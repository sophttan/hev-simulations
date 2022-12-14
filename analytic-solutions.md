# Analytical solutions for HEV simulations

## General approach
Find the theoretical probability of any secondary infection among susceptible household members $j=1...S_{i,t}$ in a household with $n$ total members given that there is 1 reported infection in household $i$ at time $t$.  

$$P(\text{any secondary infection in household i}) = 1 - P(\text{no secondary infection in household } i)$$

$$P(\text{no secondary infection in household i}) = \prod_{j=1}^{S_{i,t}} P(\text{no secondary infection in household member } j)$$


### Person-person transmission
Individuals can be infected from the community ( $I_{c,t}$ = number of current infections in the community):  

$$risk_{c} = \frac{\beta_{pp}I_{c,t}}{(1000-n)}$$

or from household contacts ( $rr$ = household relative risk):   

$$risk_{h} = \frac{rr\beta_{pp}I_{i,t}}{(n-1)}$$ 

$$P(\text{no secondary infection in household member }j) = (1-(risk_{h}+risk_{c}-risk_{h}risk_{c}))^{7}$$  


### Environmental transmission
Individuals can only be infected through contact with a contaminated environmental source:

$$risk_{e} = \beta_{e}$$  

$$P(\text{no secondary infection in household member }j) = (1-risk_{e})^{7}$$ 


### Combined person-person and environmental transmission
Individuals can be infected through person-person contact or through environmental transmission:

$$\text{risk of infection} = risk_{h}+risk_{c}+risk_{e}-risk_{h}risk_{c}-risk_{c}risk_{e}-risk_{h}risk_{e}+risk_{h}risk_{c}risk_{e}$$

$$P(\text{no secondary infection in household member }j) = (1-\text{risk of infection})^{7}$$ 
