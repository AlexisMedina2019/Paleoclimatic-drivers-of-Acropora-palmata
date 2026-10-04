## ages_kernel_reg_updt2.xlsx: same first-sheet structure as ages_kernel_reg_updt.xlsx (read by the qmd),
## with ages, groups and record inclusion taken from the audited master log (coral_ages_audit_M3.xlsx).
import pandas as pd, numpy as np
OLD='/mnt/user-data/uploads/Paleoclimatic-drivers-of-Acropora-palmata/ages_kernel_reg_updt.xlsx'
o=pd.read_excel(OLD); W='coral_ages_audit_M3.xlsx'
M=pd.read_csv('an/master_ages_M3.csv')
# --- old ID for each master-log row
pm=pd.read_excel(W,sheet_name='PM_ages'); pm=pm[pd.to_numeric(pm['ID'],errors='coerce').notna()]
sc=pd.read_excel(W,sheet_name='SC_Hubbard'); sc=sc[pd.to_numeric(sc['ID'],errors='coerce').notna()]
bz=pd.read_excel(W,sheet_name='BZ_Gischler'); bz=bz[pd.to_numeric(bz['ID'],errors='coerce').notna()]
ids=[]
for site,sh in [('Punta Maroma',pm),('Saint Croix',sc),('Belize',bz)]:
    mm=M[M.site==site]; assert len(mm)==len(sh),(site,len(mm),len(sh)); ids+=list(zip(mm.index,sh['ID'].astype(int)))
# Florida: name (spaces/asterisks removed) + stratum (= old 'source'), ties broken by nearest published age
nrm=lambda x: str(x).replace(' ','').replace('*','').strip()
fo=o[o.spot=='Florida'].copy(); fo['k']=fo.sample_site.map(nrm)
fm=M[M.site=='Florida'].copy(); fm['k']=fm.sample_id.map(nrm)
cand=fm.reset_index().merge(fo[['ID','k','source','t_ka']],left_on=['k','stratum'],right_on=['k','source'])
cand['d']=(cand.t_ka-cand.t_pub_ka).abs(); cand=cand.sort_values('d')
used_m=set(); used_o=set()
for _,r in cand.iterrows():
    if r['index'] in used_m or r.ID in used_o: continue
    ids.append((r['index'],int(r.ID))); used_m.add(r['index']); used_o.add(r.ID)
M['ID']=pd.Series(dict(ids)).reindex(M.index)
print('master rows without old ID:', M.ID.isna().sum()); print(M[M.ID.isna()][['site','sample_id','keep']].to_string())
unm=fo[~fo.ID.isin(used_o)]; print('old Florida rows not paired (duplicates):'); print(unm[['ID','sample_site','source','t_ka']].to_string())
# Florida old rows that pair with the same audited sample but were not used: they are duplicate entries -> excluded
# --- sheet 1
oi=o.set_index('ID')
rows=[]
for _,r in o.iterrows():
    if r.spot not in ('Punta Maroma','Saint Croix','Belize','Florida'):
        x=r.to_dict(); x['_status']='not audited (site not modelled)'; rows.append(x); continue
    q=M[M.ID==r.ID]
    if not len(q): rows.append({**r.to_dict(),'_status':'excluded: duplicate entry of an audited sample'}); continue
    q=q.iloc[0]
    if q.keep!='yes': rows.append({**r.to_dict(),'_status':'excluded: duplicate entry of an audited sample'}); continue
    x=r.to_dict(); x['group']=q.group; x['t_ka']=round(float(q.t_used_ka),5); x['_status']='audited'; rows.append(x)
A=pd.DataFrame(rows)
S1=A[A._status!='excluded: duplicate entry of an audited sample'][o.columns.tolist()].reset_index(drop=True)
# --- traceability sheet
T=A[['ID','sample_site','spot','group','t_ka','_status']].rename(columns={'group':'group_updt2','t_ka':'t_ka_updt2','_status':'status'})
T=T.merge(o[['ID','group','t_ka']].rename(columns={'group':'group_updt1','t_ka':'t_ka_updt1'}),on='ID')
T.loc[T.status.str.startswith('excluded'),['t_ka_updt2','group_updt2']]=np.nan
ML=M.drop(columns=['site']).rename(columns={'sample_id':'sample_id_audit'})
T=T.merge(ML,on='ID',how='left')
T['dt_ka']=(T.t_ka_updt2-T.t_ka_updt1).round(5)
T['bin_updt1']=np.maximum(0.1,np.ceil(np.round(T.t_ka_updt1,6)*10)/10)
T['bin_updt2']=np.where(T.t_ka_updt2.notna(),np.maximum(0.1,np.ceil(np.round(T.t_ka_updt2.fillna(0),6)*10)/10),np.nan)
T['flags']=T['flags'].astype(object).where(T['flags'].notna(),'')
T=T[['ID','spot','sample_site','sample_id_audit','sample_id_source','status','group_updt1','group_updt2','t_ka_updt1','t_ka_updt2','dt_ka','bin_updt1','bin_updt2','method','source','stratum','subregion','t_pub_ka','ref_year','err_2s_ka','conv14C','conv14C_err','dR','dR_err','curve','t_LT_marine20_ka','t_sens_dR_regional_ka','keep','in_model','flags']]
T.to_csv('updt2_trace.csv',index=False); S1.to_csv('updt2_registros.csv',index=False)
pal=lambda d,c: d[(d.group=='palmata')&(d[c]>-0.1)&(d[c]<=6)]
print('sheet1 rows', len(S1)); print(S1.groupby(['spot','group']).size())
print('palmata 0-6 by site (updt2):'); print(pal(S1,'t_ka').groupby('spot').size())
print('in_model master:', (M.in_model=='yes').sum())
