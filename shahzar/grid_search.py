import numpy as np
import pandas as pd
from scipy import stats

eps = np.finfo(float).eps

time = 365
num_weeks_inc = 28
num_weeks_inf = 7

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

def SEIR(beta_H, beta_C, inc, inf):
    hh_size = create_hh()
    
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
    
    results = data.loc[:, 'ID':'HH']
    results['TYPE'] = np.nan
    results['TIME'] = np.nan
    
    for t in range(time):
        if t % 10 == 0: print(t, end = ' ')
        recovered = (data['INF'] > 0) & (data['I_count'] == data['INF'])
        if sum(recovered) > 0:
            data.loc[recovered, 'R'] = 1
            data.loc[recovered, 'I'] = 0

            # Why is this count set to 0 once they recover?
            data.loc[recovered, 'I_count'] = 0

        new_inf = (data['INC'] > 0) & (data['E_count'] == data['INC'])
        num_new_inf = sum(new_inf)

        if num_new_inf > 0:
            # Why do new infections have SD = 1?
            random_inf = np.round(stats.norm.rvs(inf, 1, size = num_new_inf))
            data.loc[new_inf, 'I'] = 1
            data.loc[new_inf, 'INF'] = random_inf
            data.loc[new_inf, 'E'] = 0

            # Why is this set to 0?
            data.loc[new_inf, 'E_count'] = 0

        I_H = data.groupby('HH').sum()['I']
        summary = pd.DataFrame({'I_H': I_H, 
                                'I_C': data.sum()['I'] - I_H
                               })
        dd = data[['HH', 'S']].copy()

        dd['I_H'] = dd.apply(lambda x: summary.loc[x['HH'], 'I_H'], axis = 1)
        dd['I_C'] = dd.apply(lambda x: summary.loc[x['HH'], 'I_C'], axis = 1)

        risk_H = dd['S'] * beta_H * dd['I_H'] / 1000
        risk_C = dd['S'] * beta_C * dd['I_C'] / 1000

        new_inf_H = stats.binom.rvs(1, risk_H, size = data.shape[0])
        new_inf_C = stats.binom.rvs(1, risk_C, size = data.shape[0])
        new_exposed = (new_inf_H == 1) | (new_inf_C == 1)

        num_new_exposed = sum(new_exposed)
        if num_new_exposed > 0:
            data.loc[new_exposed, 'E'] = 1
            # Why is this calculated only upon exposure?
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
    incidence = state.shape[0]/1000
    
    prop_hh = 0
    if state.shape[0] != 0:
        prop_hh = state[state['TYPE'] == 'H'].shape[0] / state.shape[0]
    
    return incidence, prop_hh

def mse(results, target):
    return ((np.array(metrics(results)) - target)**2).sum()

def bce(results, target):
    incidence, prop_hh = metrics(results)
    
    incidence_loss = target[0] * np.log(incidence + eps) + (1 - target[0]) * np.log(1 - incidence + eps)
    prop_hh_loss = target[1] * np.log(prop_hh + eps) + (1 - target[1]) * np.log(1 - prop_hh + eps)
    return incidence_loss + prop_hh_loss

vals = np.arange(0.05, 1, 0.05)
target = np.array([0.1, 0.25])

n = len(vals)
bces = np.zeros((n, n))
mses = np.zeros((n, n))
for i in range(n):
    beta_H = vals[i]
    for j in range(n):
        beta_C = vals[j]
        results = SEIR(beta_H, beta_C, num_weeks_inc, num_weeks_inf)
        
        bces[i][j] = bce(results, target)
        mses[i][j] = mse(results, target)

np.save('bces.npy', bces)
np.save('mses.npy', mses)
np.save('vals.npy', vals)