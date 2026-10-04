## Age draws for the M3 Monte Carlo. One row per in-model record, N columns.
## 14C: sampled from the Marine20 calibrated PDF (conv14C, conv14C_err, dR, dR_err), then shifted so that its
##      median equals t_used_ka (the published / audited central value; offsets are <=1 yr except Looe Key and a
##      few old Florida samples, 18-52 yr, where LT's calibration software differs slightly).
## U-Th: normal(t_used_ka, err_2s_ka/2).
## Belize sensitivity set: 14C with regional dR (-114 +/- 70) centred on t_sens_dR_regional_ka.
import numpy as np, pandas as pd, sys
sys.path.insert(0, '..'); from fastcal import C
N = 1000; rng = np.random.default_rng(20261002)
M = pd.read_csv('master_ages_M3.csv'); m = M[M.in_model == 'yes'].reset_index(drop=True)
fl_dre = m.loc[(m.site == 'Florida') & (m.method == '14C'), 'dR_err'].median()
miss = m.dR_err.isna() & (m.method == '14C'); log = m.loc[miss, 'sample_id'].tolist()
m.loc[miss, 'dR_err'] = fl_dre
g, mu, s = C['marine20']
def pdf_draw(c, e, dr, dre, centre):
    sd2 = e**2 + dre**2 + s**2; p = np.exp(-(c - dr - mu)**2 / (2*sd2)) / np.sqrt(sd2); p /= p.sum()
    cd = np.cumsum(p); med = g[np.searchsorted(cd, .5)]
    x = rng.choice(g, size=N, p=p) + rng.uniform(-0.5, 0.5, N)
    return (x - med) / 1000 + centre, (med/1000 - centre)*1000
D = np.zeros((len(m), N)); DB = np.zeros((len(m), N)); off = []
for i, r in m.iterrows():
    if r.method == '14C':
        D[i], o = pdf_draw(r.conv14C, r.conv14C_err, r.dR, r.dR_err, r.t_used_ka); off.append((r.site, r.sample_id, o))
        if r.site == 'Belize':
            DB[i], _ = pdf_draw(r.conv14C, r.conv14C_err, -114, 70, r.t_sens_dR_regional_ka)
        else: DB[i] = D[i]
    else:
        D[i] = rng.normal(r.t_used_ka, max(r.err_2s_ka, 1e-6)/2, N); DB[i] = D[i]
pd.DataFrame(D).round(5).to_csv('draws_ages.csv', index=False)
pd.DataFrame(DB).round(5).to_csv('draws_ages_BZdR.csv', index=False)
m[['site','sample_id','stratum','method','t_used_ka','err_2s_ka','t_sens_dR_regional_ka']].to_csv('draws_meta.csv', index=False)
o = pd.DataFrame(off, columns=['site','sample_id','offset_yr']); o.to_csv('draws_offsets.csv', index=False)
print('dR_err filled with Florida median', round(fl_dre,1), 'for', log)
print(o.groupby('site').offset_yr.describe().round(1))
print('draw SD (yr) by site/method:'); m['sd'] = D.std(1)*1000; print(m.groupby(['site','method']).sd.describe()[['mean','min','max']].round(0))
