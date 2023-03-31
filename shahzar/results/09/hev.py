import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
import pandas as pd
from scipy import stats
import timeit

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
    results.loc[0, 'TYPE'] = 'I' # index case is Type 0 
    results['TIME'] = np.nan
    #results['TIME_I'] = np.nan
    #results['TIME_E'] = np.nan
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
            data.loc[new_exposed, 'E'] = 1
            data.loc[new_exposed, 'INC'] = random_inc
            
            #print(np.sum( (new_inf_H == 1) & (~pd.isna(results['TYPE'])) ))
            # The NA checks may not be necessary, since S = 0 for anyone who already has been assigned a TYPE
            results['TYPE'].where(~(new_inf_H == 1) | ~pd.isna(results['TYPE']), 'H', inplace = True)
            results['TYPE'].where(~(new_inf_C == 1) | ~pd.isna(results['TYPE']), 'C', inplace = True)
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
        data.loc[data['E'] == 1, 'S'] = 0

    return results

def metrics(results):
    state = results[~pd.isna(results['TIME'])]
    idc = state.shape[0]/1000
    
    sar = np.nan
    if idc != 0:
        sar = np.mean(results['I_num'] / results['S_num'])
    
    return idc, sar

def score(results, target):
    return np.sum((np.array(metrics(results)) - target)**2)

beta_Hs = np.linspace(5, 35, 61)
beta_Cs = np.linspace(0, 1, 21)

cc, hh = np.meshgrid(beta_Cs, beta_Hs)
a, b = cc.shape

np.save('cc.npy', cc)
np.save('hh.npy', hh)

reps = 25
sars = np.zeros((a, b, reps))
idcs = np.zeros((a, b, reps))
tocs = np.zeros(a * b * reps)
for i in range(a):
    for j in range(b):
        beta_H = hh[i][j]
        beta_C = cc[i][j]
        
        idc_list = np.array([None] * reps)
        sar_list = np.array([None] * reps)
        for k in range(reps):
            t_0 = timeit.default_timer()
            print('{}/{}'.format(beta_H, 35), '\t', '{}/{}'.format(beta_C, 1), '\t', '{}/{}'.format(k, reps), end = '\t')
            results = SEIR(beta_H, beta_C, inc, inf, verbose = 0)
            idc, sar = metrics(results)
            idc_list[k] = idc
            sar_list[k] = sar
            print(np.round(idc, 3), '\t', np.round(sar, 3), end = '\t')
            t_1 = timeit.default_timer()
            print(round((t_1 - t_0), 3))
            tocs[b*reps*i + reps*j + k] = t_1 - t_0

        idcs[i][j] = idc_list
        sars[i][j] = sar_list
        
        np.save('idcs.npy', idcs)
        np.save('sars.npy', sars)
        np.save('tocs.npy', tocs)
        print(idc_list)
        print(sar_list)
        print()

np.save('idcs.npy', idcs)
np.save('sars.npy', sars)
np.save('tocs.npy', tocs)