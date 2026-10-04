## Gap robustness: for each site, the record-free intervals of the corrected central ages that overlap the reported gaps;
## under the age Monte Carlo, the length of the interval between the bracketing records and whether any record intrudes.
import pandas as pd, numpy as np
m=pd.read_csv('draws_meta.csv'); D=pd.read_csv('draws_ages.csv').values
G={'Punta Maroma':[(0.71,1.44)],'Belize':[(1.99,2.72),(3.68,4.20),(5.54,6.0)]}
out=[]
for s,gs in G.items():
    idx=np.where(m.site==s)[0]; t=m.t_used_ka.values[idx]; o=np.argsort(t); ts=t[o]
    gaps=np.diff(ts)
    for a,b in gs:
        # largest record-free interval overlapping the reported gap
        cand=[(i,gaps[i]) for i in range(len(gaps)) if ts[i]<b and ts[i+1]>a]
        i,_=max(cand,key=lambda x:x[1]); y,old=idx[o[i]],idx[o[i+1]]
        d=D[idx]; dy=D[y]; do=D[old]
        L=(do-dy); intr=((d>np.minimum(dy,do)[None,:])&(d<np.maximum(dy,do)[None,:])).sum(0)-0
        out.append(dict(site=s,reported=f'{a}-{b}',central_gap=f'{ts[i]:.3f}-{ts[i+1]:.3f}',central_len_ka=round(ts[i+1]-ts[i],3),
            young=m.sample_id[y]+' '+m.method[y],old=m.sample_id[old]+' '+m.method[old],
            len_q025=round(np.quantile(L,.025),3),len_med=round(np.median(L),3),len_q975=round(np.quantile(L,.975),3),
            p_no_intrusion=round((intr==0).mean(),3),p_len_gt_0_3=round((L>0.3).mean(),3)))
r=pd.DataFrame(out); r.to_csv('gaps_mc.csv',index=False); print(r.to_string())
