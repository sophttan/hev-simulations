import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
import pandas as pd
from scipy import stats

time = 365
inc = 28
inf = 7

def create_hh():
    hh_size = np.random.choice([3, 4, 5, 6], size = 340, replace = True)
    hh_size = hh_size[np.cumsum(hh_size) < 1000]
    
    leftover = 1000 - hh_size.sum()
    if leftover < 3:
        # Randomly sample leftover amount indices from hh_size, and add one
        hh = np.arange(len(hh_size))
        idx = np.random.choice(hh[hh_size < 6], size = leftover, replace = False)
        hh_size[idx] += 1
    else:
        hh_size = np.append(hh_size, leftover)
    return hh_size

def SEIR(beta_H, beta_C, inc, inf, verbose = 0):
    hh_size = create_hh()
    
    # ID: ID of individual
    # SIZE: size of individual's household
    # HH: ID of individual's household
    # S: susceptibility status
    # E: exposed status
    # E_count: number of days since exposed
    # I: infectious status
    # I_count: number of days since infectious
    # R: recovered status
    # INC: incubation period
    # INF: infectious period
    # S_num: number of susceptible people in household when infectious period begins
    # I_num: number of people that this person infected over the infectious period
    data = pd.DataFrame({'ID': range(1000),
                         'SIZE': np.repeat(hh_size, repeats = hh_size),
                         'HH': np.repeat(range(len(hh_size)), repeats = hh_size),
                         'S': np.append(0, np.ones(999)),
                         'E': np.append(1, np.zeros(999)),
                         'E_count': np.append(1, np.zeros(999)),
                         'I': np.zeros(1000),
                         'I_count': np.zeros(1000),
                         'R': np.zeros(1000),
                         'INC': np.append(np.round(stats.norm.rvs(inc, 2)), np.zeros(999)),
                         'INF': np.zeros(1000)
                        })
    
    # Create frame for storing results
    results = data.loc[:, 'ID':'HH']
    results['TYPE'] = np.nan
    results.loc[0, 'TYPE'] = '0' # index case is Type 0 
    results['TIME'] = np.nan
    results['S_num'] = np.nan
    results['I_num'] = 0
    
    for t in range(time):
        if verbose:
            if t % 10 == 0: print(t, end = ' ')
        
        # Anyone who has been infectious for as many days as their infectious
        # period is now recovered.
        recovered = (data['INF'] > 0) & (data['I_count'] == data['INF'])
        if sum(recovered) > 0:
            data.loc[recovered, 'R'] = 1
            data.loc[recovered, 'I'] = 0
            data.loc[recovered, 'I_count'] = 0

        # Anyone who has been incubating for as many days as their incubation
        # period is now infectious.
        new_inf = (data['INC'] > 0) & (data['E_count'] == data['INC'])
        num_new_inf = sum(new_inf)
        if num_new_inf > 0:
            random_inf = np.round(stats.norm.rvs(inf, 1, size = num_new_inf))
            data.loc[new_inf, 'I'] = 1
            data.loc[new_inf, 'INF'] = random_inf
            data.loc[new_inf, 'E'] = 0
            data.loc[new_inf, 'E_count'] = 0
            
            # If newly infected and there isn't already a time recorded, then update with the time t.
            results['TIME'].where(~(new_inf == 1) | ~pd.isna(results['TIME']), t, inplace = True)
            #results['TIME_I'].where(~(new_inf == 1) | ~pd.isna(results['TIME_I']), t, inplace = True)
            
            # Get the number of susceptible people in each household
            S_num = data.groupby('HH').sum()

            # Get the households of the newly infectious people
            hh_idx = data.loc[new_inf, :]['HH']
            
            # Save in results the number of susceptible people
            results.loc[new_inf, 'S_num'] = S_num.loc[hh_idx, 'S'].values
            
        # I_H is the number of infections in each household.
        # I_C is the number of infections outside a given household.
        I_H = data.groupby('HH').sum()['I']
        summary = pd.DataFrame({'I_H': I_H, 
                                'I_C': data.sum()['I'] - I_H
                               })
        # dd is a frame where each individual is assigned their household's
        # I_H and I_C numbers.
        dd = data[['HH', 'S']].copy()
        dd['I_H'] = dd.apply(lambda x: summary.loc[x['HH'], 'I_H'], axis = 1)
        dd['I_C'] = dd.apply(lambda x: summary.loc[x['HH'], 'I_C'], axis = 1)

        # Calculate household risk and community risk
        risk_H = dd['S'] * beta_H * dd['I_H'] / 1000
        risk_C = dd['S'] * beta_C * dd['I_C'] / 1000

        # Calculate new household and community infections
        new_inf_H = stats.binom.rvs(1, risk_H, size = data.shape[0])
        new_inf_C = stats.binom.rvs(1, risk_C, size = data.shape[0])
        new_exposed = (new_inf_H == 1) | (new_inf_C == 1)

        num_new_exposed = sum(new_exposed)
        if num_new_exposed > 0:
            data.loc[new_exposed, 'E'] = 1
            random_inc = np.round(stats.norm.rvs(inc, 2, size = num_new_exposed))
            data.loc[new_exposed, 'INC'] = random_inc
            data.loc[new_exposed, 'S'] = 0
            
            #print(np.sum( (new_inf_H == 1) & (~pd.isna(results['TYPE'])) ))
            # The NA checks may not be necessary, since S = 0 for anyone who already has been assigned a TYPE
            results['TYPE'].where(~(new_inf_C == 1) | ~pd.isna(results['TYPE']), 'C', inplace = True)
            results['TYPE'].where(~(new_inf_H == 1) | ~pd.isna(results['TYPE']), 'H', inplace = True)
            #results['TIME_E'].where(~(new_exposed == 1) | ~pd.isna(results['TIME_E']), t, inplace = True)
            
            # Number of new infections in each household
            rr = results.loc[new_inf_H == 1, :].groupby('HH')['TYPE'].count()
            
            # Get individuals with the smallest infectious counts
            ids = data[data['I'] == 1].groupby('HH').min()['ID']
            
            # Increase I_num by the amount of people in household that got a H-type infection
            res_idx = results['ID'].isin(ids) & (results['HH'].isin(rr.index))
            results.loc[res_idx, 'I_num'] += rr.values
            
            # After the above analysis, assign 'B' to people who were infected by both.
            results['TYPE'].where(~((new_inf_H == 1) & (new_inf_C == 1)) | ~pd.isna(results['TYPE']), 'B', inplace = True)
            
        data.loc[data['E'] == 1, 'E_count'] += 1
        data.loc[data['I'] == 1, 'I_count'] += 1

    return results

def metrics(results):
    idc = np.mean(~pd.isna(results['TYPE']))
    
    sar = np.nan
    if idc != 0:
        sar = np.mean(results['I_num'] / results['S_num'])
    
    return idc, sar

def score(results, target):
    return np.sum((np.array(metrics(results)) - target)**2)

#### Metropolis algorithm ####

# Create likelihood from the score of the state
def likelihood(state):
    beta_H, beta_C = state
    
    liks = [None] * N
    for i in range(N):
        results = SEIR(beta_H, beta_C, inc, inf)
        liks[i] = -np.log(score(results, target))
    return np.mean(liks)

# Proposal function
def q(state):
    beta_H, beta_C = state
    r = [beta_H**2, beta_C**2 / 0.0001]
    v = np.array([beta_H, beta_C / 0.0001])
    return np.maximum(stats.gamma.rvs(r, scale = 1 / v, size = 2), 1e-3)

# MCMC
def metropolis(start, num_iter):
    chain = np.zeros((num_iter + 1, 2))
    liks = np.zeros((num_iter + 1, 2))
    
    # Initialize current state.
    curr = start
    curr_lik = likelihood(curr)
    
    # Initialize best state.
    best = curr
    best_lik = curr_lik
    for i in range(num_iter):
        # Save the current state and its likelihood.
        chain[i] = curr
        liks[i] = curr_lik
        
        # Print current state and likelihood.
        print(i, '\t', curr, end = ' ')
        print("%0.3f" % curr_lik, end = '\t')
        
        # Get a proposed state and calculate its likelihood.
        prop = q(curr)
        print(prop, end = ' ')
        
        prop_lik = likelihood(prop)
        print("%0.3f" % prop_lik, end = '\t')
        
        # Compute the ratio of the scores of the two states
        # and flip a coin.
        r = np.exp(prop_lik - curr_lik)
        p = stats.uniform.rvs()
        print("%0.3f" % r, ' ', "%0.3f" % p)
        # Transition if the proposed state is better or
        # if the coin flip succeeds.
        if p < r:
            curr = prop
            curr_lik = prop_lik
            
            # If the new likelihood is better than the
            # best we've seen so far, replace the best.
            if curr_lik > best_lik:
                best = curr
                best_lik = curr_lik
        
        # Save the chain, best state, and likelihoods
        # so far.
        np.save('chain.npy', chain)
        np.save('liks.npy', liks)
        np.save('best.npy', best)
    
    chain[num_iter] = curr
    liks[num_iter] = curr_lik

    np.save('chain.npy', chain)
    np.save('liks.npy', liks)
    np.save('best.npy', best)
    return chain, liks, best

# Solve for optimal values via MCMC
target = np.array([0.3, 0.25]) # target values
N = 5 # number of times over which to average likelihood

#chain, liks, best = metropolis([30, 0.12], 1000)
#np.save('chain.npy', chain)
#np.save('liks.npy', liks)
#np.save('best.npy', best)