import pandas as pd, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter as L
F='Arial'; NORM=Font(name=F); BOLD=Font(name=F,bold=True); HDR=PatternFill('solid',fgColor='D9E1F2'); RED=PatternFill('solid',fgColor='F8CBAD'); GRY=PatternFill('solid',fgColor='EDEDED')
S1=pd.read_csv('updt2_registros.csv'); T=pd.read_csv('updt2_trace.csv')
T.loc[T.status.str.startswith('excluded') & T.sample_id_audit.isna(),'flags']='entrada repetida en updt1 (tercera copia FKRT de una muestra DRTO, Stathakopoulos et al. 2025): se cuenta una vez'
wb=Workbook(); ws=wb.active; ws.title='Registros'
def put(ws,df,widths=None):
    for j,c in enumerate(df.columns,1): x=ws.cell(1,j,c); x.font=BOLD; x.fill=HDR; x.alignment=Alignment(wrap_text=True)
    for i,row in enumerate(df.itertuples(index=False),2):
        for j,v in enumerate(row,1):
            if isinstance(v,float) and np.isnan(v): v=None
            ws.cell(i,j,v).font=NORM
    ws.freeze_panes='B2'
    for k,v in (widths or {}).items(): ws.column_dimensions[L(k)].width=v
put(ws,S1,{1:6,2:18,3:34,4:13,5:22,6:9,7:13,8:9})
ws2=wb.create_sheet('Trazabilidad'); put(ws2,T,{1:6,2:13,3:18,4:18,5:18,6:30,7:12,8:12,9:9,10:9,11:8,12:7,13:7,14:7,15:28,16:30,17:18,18:8,19:7,20:8,21:8,22:7,23:8,24:6,25:9,26:9,27:9,28:5,29:7,30:90})
for i,s in enumerate(T.status,2):
    if s.startswith('excluded'):
        for j in range(1,T.shape[1]+1): ws2.cell(i,j).fill=RED
    elif s.startswith('not audited'):
        for j in range(1,T.shape[1]+1): ws2.cell(i,j).fill=GRY
ws3=wb.create_sheet('README')
txt=[('ages_kernel_reg_updt2.xlsx — cronología de corales usada en regional_analysis_Final.qmd',BOLD),
('Creado 2026-10-03. Sustituye a ages_kernel_reg_updt.xlsx (updt1) como entrada del qmd (chunk datos-reg).',NORM),('',NORM),
('Hoja Registros (la que lee el qmd): mismas 8 columnas y mismos ID que updt1. Cambian:',BOLD),
(' • t_ka = edad auditada en ka BP (1950): U-Th referidas a 1950; ¹⁴C recalibradas con Marine20 y ΔR local (Florida: valores de L. T. Toth).',NORM),
(' • group: 7 A. cervicornis de Belice codificadas antes como palmata se reclasifican (Other species).',NORM),
(' • Se eliminan 22 filas que eran entradas repetidas de una misma muestra (Florida 17, Saint Croix 5).',NORM),
(' • Una edad posterior a 1950 (Punta Maroma SURF7(a), −0.00568 ka) se conserva con su valor; el qmd la asigna al bin 0–0.1 ka.',NORM),
(' • Barbados (105 filas) no se modela y se copia sin cambios (no auditado).',NORM),('',NORM),
('Hoja Trazabilidad: una fila por registro de updt1 (enlazada por ID) con su estado (audited / excluded / not audited), la edad y el grupo en updt1 y updt2, el bin en cada versión, el método, la fuente original, la edad publicada, el año de referencia, el error 2σ, la edad convencional ± σ, el ΔR ± σ, la curva, la edad de LT (Florida), la edad con ΔR regional (Belice) y las notas de auditoría.',NORM),('',NORM),
('Fuente de todas las correcciones: coral_ages_audit_M3.xlsx (hojas por sitio, Master_log, Changelog) y master_ages_M3.csv. Construido con build_updt2.py y write_updt2.py.',NORM),('',NORM),
('Resumen A. palmata 0–6 ka (filtro del qmd: t_ka > −0.1 y ≤ 6):',BOLD)]
pal=S1[(S1.group=='palmata')&(S1.t_ka>-0.1)&(S1.t_ka<=6)]
o1=T[(T.group_updt1=='palmata')&(T.t_ka_updt1>=0)&(T.t_ka_updt1<=6)]
for s in ['Punta Maroma','Florida','Belize','Saint Croix']:
    txt.append((f' • {s}: updt1 {int((o1.spot==s).sum())} → updt2 {int((pal.spot==s).sum())}',NORM))
txt.append((f' • Total: updt1 {len(o1[o1.spot!="Barbados"])} → updt2 {len(pal[pal.spot!="Barbados"])}',NORM))
for i,(t,f) in enumerate(txt,1): c=ws3.cell(i,1,t); c.font=f
ws3.column_dimensions['A'].width=160
wb.active=0
wb.save('ages_kernel_reg_updt2.xlsx')
print([ (t) for t,_ in txt[-5:]])
