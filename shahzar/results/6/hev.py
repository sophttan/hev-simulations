import numpy as np
import pandas as pd
import timeit
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
                         'INF': np.zeros(1000),
                        })
    
    # Create frame for storing results
    results = data.loc[:, 'ID':'HH']
    results['TYPE'] = np.nan
    results.loc[0, 'TYPE'] = '0'
    results['TIME'] = np.nan
    results.loc[0, 'TIME'] = 0
    
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
            data.loc[new_exposed, 'E'] = 1
            data.loc[new_exposed, 'INC'] = random_inc

            results['TYPE'].where(~((new_inf_H == 1) & (new_inf_C == 1)) | ~pd.isna(results['TYPE']), 'B', inplace = True)
            results['TYPE'].where(~(new_inf_H == 1) | ~pd.isna(results['TYPE']), 'H', inplace = True)
            results['TYPE'].where(~(new_inf_C == 1) | ~pd.isna(results['TYPE']), 'C', inplace = True)
            results['TIME'].where(~(new_exposed == 1) | ~pd.isna(results['TIME']), t, inplace = True)
        
        data.loc[data['E'] == 1, 'E_count'] += 1
        data.loc[data['I'] == 1, 'I_count'] += 1
        data.loc[data['E'] == 1, 'S'] = 0

    return results

def metrics(results):
    state = results[~pd.isna(results['TIME'])]
    idc = state.shape[0]/1000
    
    #p_hh = np.nan
    #if i != 0:
    #    p_hh = state[state['TYPE'] == 'H'].shape[0] / state.shape[0]
    #return i, p_hh
    
    sar = np.nan
    if idc != 0:
        num_primary = np.sum(results.groupby('HH')['TIME'].sum() > 0) # households that were infected
        idx = results.groupby('HH')['TYPE'].apply(lambda x: ~np.all(x.isna()))
        num_contact = (results.groupby('HH')['SIZE'].sum()[idx]**(1/2)).sum() # total people in all those households
        sar = state[(state['TYPE'] == 'H') | (state['TYPE'] == 'B')].shape[0] / (num_contact - num_primary)
    
    return idc, sar

def score(results, target):
    return np.sum((np.array(metrics(results)) - target)**2)

beta_Hs = np.linspace(10, 50, 41)
beta_Cs = np.linspace(0, 2, 41)

cc, hh = np.meshgrid(beta_Cs, beta_Hs)
a, b = cc.shape

np.save('cc_5.npy', cc)
np.save('hh_5.npy', hh)

reps = 1
sars = np.zeros((a, b))
idcs = np.zeros((a, b))
for i in range(a):
    for j in range(b):
        beta_H = hh[i][j]
        beta_C = cc[i][j]
        
        t_0 = timeit.default_timer()
        print('{}/{}'.format(i, a), j, end = '\t')
        results = SEIR(beta_H, beta_C, inc, inf, verbose = 1)
        idc, sar = metrics(results)
        t_1 = timeit.default_timer()
        print(round((t_1 - t_0), 3), end = ' ')
        idcs[i][j] = idc
        sars[i][j] = sar
        
        np.save('idcs_5.npy', idcs)
        np.save('sars_5.npy', sars)
        print()

np.save('idcs_5.npy', idcs)
np.save('sars_5.npy', sars)