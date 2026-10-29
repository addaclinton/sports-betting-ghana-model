"""Reproduce all numerical results and figure data of the manuscript.

Usage:  python scripts/run_analysis.py [--outdir DIR]
Writes fig2.csv, fig3.csv, fig4.csv and summary.json (default: ./data).
Runtime: about a minute. Latin hypercube seed is fixed (20261027).
"""
import argparse, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
_ap = argparse.ArgumentParser()
_ap.add_argument('--outdir', default=str(ROOT / 'data'))
OUT = Path(_ap.parse_args().outdir); OUT.mkdir(parents=True, exist_ok=True)

import numpy as np, csv, json
from model import *
from scipy.optimize import brentq
from scipy.stats import qmc, rankdata

names=['beta','a1','a2','g1','g2','g3','mu']
def R0_of(p): return R0f(p['beta'],p['a1'],p['a2'],p['g1'],p['g2'],p['g3'],p['mu'])
base=R0_of(P)
sens={}
for n in names:
    h=1e-6*P[n]; pp=dict(P); pm=dict(P); pp[n]+=h; pm[n]-=h
    sens[n]=(R0_of(pp)-R0_of(pm))/(2*h)*P[n]/base

t,Y=rk4(y0,P,T)
with open(OUT/'fig2.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t','S','E','R','P','Q'])
    for ti,yi in zip(t,Y): w.writerow([f"{ti:.2f}"]+[f"{v/1e6:.6f}" for v in yi])
sc={'base':P,'beta30':{**P,'beta':P['beta']*0.7},'g3_50':{**P,'g3':P['g3']*1.5}}
res={k:rk4(y0,p,T)[1] for k,p in sc.items()}
with open(OUT/'fig3.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['t','base','beta30','g3_50'])
    for i,ti in enumerate(t): w.writerow([f"{ti:.2f}"]+[f"{res[k][i,3]/1e3:.4f}" for k in ['base','beta30','g3_50']])
b30=res['beta30'][:,3]; pk=int(np.argmax(b30))

def eqp(R0t):
    p=dict(P); p['beta']=P['beta']*R0t/base; _,e=endemic(p); return e
with open(OUT/'fig4.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['R0','prev','prob'])
    for r in np.linspace(1.0,6.0,101):
        e=eqp(r); w.writerow([f"{r:.3f}",f"{e[1:4].sum()/e.sum():.5f}",f"{e[3]/e.sum():.5f}"])
R0_tgm=brentq(lambda r:eqp(r)[1:4].sum()/eqp(r).sum()-TGM_PREV,1.01,50)
e_t=eqp(R0_tgm)
Rs,eq=endemic(P)

# TGM-matched scenario, simulated
Ptgm={**P,'beta':P['beta']*R0_tgm/base}
tg={'base':Ptgm,'beta30':{**Ptgm,'beta':Ptgm['beta']*0.7},'g3_50':{**Ptgm,'g3':Ptgm['g3']*1.5}}
tgres={k:rk4(y0,p,T)[1] for k,p in tg.items()}

# ---------- uncertainty analysis (vectorised RK4, fixed bug: a1,a2 held at baseline in scenarios) ----------
rng_lo=dict(R0=0.8,g1=0.125,g2=0.025,g3=0.025,delta=0.034)
rng_hi=dict(R0=6.0,g1=2.0,g2=0.4,g3=0.54,delta=0.88)
keys=['R0','g1','g2','g3','delta']
n=20000
u=qmc.LatinHypercube(d=5,seed=20261027).random(n)
X={k:np.exp(np.log(rng_lo[k])+u[:,i]*(np.log(rng_hi[k])-np.log(rng_lo[k]))) for i,k in enumerate(keys)}
def vec_pars(X):
    c=X['g3']+mu
    a2=(share['P']/share['R'])*c; b=a2+X['g2']+mu; a1=(share['R']/share['E'])*b
    a=a1+X['g1']+mu
    beta=X['R0']/(1/a+a1/(a*b)+a1*a2/(a*b*c))
    return dict(beta=beta,a1=a1,a2=a2,g1=X['g1'],g2=X['g2'],g3=X['g3'],delta=X['delta'])
def simulate(pp,T=30.0,h=0.05):
    m=len(pp['beta']); Y=np.zeros((m,5)); Y[:]=y0
    def f(Y):
        S,E,R,Pp,Q=Y.T; N=Y.sum(1); lam=pp['beta']*(E+R+Pp)/N
        return np.stack([Lam+pp['delta']*Q-lam*S-mu*S, lam*S-(pp['a1']+pp['g1']+mu)*E,
          pp['a1']*E-(pp['a2']+pp['g2']+mu)*R, pp['a2']*R-(pp['g3']+mu)*Pp,
          pp['g1']*E+pp['g2']*R+pp['g3']*Pp-(pp['delta']+mu)*Q],1)
    for _ in range(int(round(T/h))):
        k1=f(Y);k2=f(Y+h/2*k1);k3=f(Y+h/2*k2);k4=f(Y+h*k3);Y=Y+h/6*(k1+2*k2+2*k3+k4)
    return Y
# --- validation: vectorised solver at baseline values must reproduce the scalar solver ---
Xb={k:np.array([v]) for k,v in dict(R0=1.3,g1=g1,g2=g2,g3=g3,delta=delta).items()}
vb=simulate(vec_pars(Xb))[0]
check=float(np.max(np.abs(vb-res['base'][-1])/np.maximum(res['base'][-1],1)))
vp=vec_pars(Xb); assert abs(vp['beta'][0]-P['beta'])<1e-12 and abs(vp['a1'][0]-P['a1'])<1e-12

bp=vec_pars(X); P30=simulate(bp)[:,3]
pb=dict(bp); pb['beta']=bp['beta']*0.7; P30b=simulate(pb)[:,3]
pg=dict(bp); pg['g3']=bp['g3']*1.5; P30g=simulate(pg)[:,3]
mask=X['R0']>1
red_b=1-P30b/P30; red_g=1-P30g/P30
def pct(a,q): return float(np.percentile(a,q))
def prcc(Xm,yv):
    R=np.column_stack([rankdata(Xm[:,i]) for i in range(Xm.shape[1])]); ry=rankdata(yv)
    out=[]
    for i in range(R.shape[1]):
        Z=np.column_stack([np.ones(len(ry)),np.delete(R,i,1)])
        r1=R[:,i]-Z@np.linalg.lstsq(Z,R[:,i],rcond=None)[0]; r2=ry-Z@np.linalg.lstsq(Z,ry,rcond=None)[0]
        out.append(float(np.corrcoef(r1,r2)[0,1]))
    return out
Xm=np.column_stack([np.log(X[k]) for k in keys])
pr=prcc(Xm,P30)
# also: share of R0>1 draws in which beta cut brings R0 below 1 ; g3 boost does
def R0vec(pp): 
    a=pp['a1']+pp['g1']+mu; b=pp['a2']+pp['g2']+mu; c=pp['g3']+mu
    return pp['beta']*(1/a+pp['a1']/(a*b)+pp['a1']*pp['a2']/(a*b*c))
below_b=float((R0vec(pb)[mask]<1).mean()); below_g=float((R0vec(pg)[mask]<1).mean())
unc=dict(n=n,check_vs_scalar=check,share_gt1=float(mask.mean()),
  P30_med=pct(P30[mask],50),P30_lo=pct(P30[mask],5),P30_hi=pct(P30[mask],95),
  redb_med=pct(red_b[mask],50),redb_lo=pct(red_b[mask],5),redb_hi=pct(red_b[mask],95),
  redg_med=pct(red_g[mask],50),redg_lo=pct(red_g[mask],5),redg_hi=pct(red_g[mask],95),
  frac_beta_better=float((red_b[mask]>red_g[mask]).mean()),
  below1_beta=below_b,below1_g3=below_g,
  prcc=dict(zip(keys,pr)),ranges=dict(lo=rng_lo,hi=rng_hi))
out=dict(sens=sens,R0=base,beta=P['beta'],a1=P['a1'],a2=P['a2'],mu=mu,Lam=Lam,
  eq=eq.tolist(),end30=Y[-1].tolist(),
  scen={k:dict(R0=R0_of(p),P30=float(res[k][-1,3]),bet30=float(res[k][-1,1:4].sum())) for k,p in sc.items()},
  b30_peak=float(b30[pk]),b30_peak_t=float(t[pk]),
  beta_thr=P['beta']/base,g3_for_R0_1=float(brentq(lambda x:R0_of({**P,'g3':x})-1,P['g3'],50)),
  R0_tgm=float(R0_tgm),tgm_beta=float(Ptgm['beta']),tgm_eq_share=(e_t/e_t.sum()).tolist(),
  tgm_scen={k:dict(R0=R0_of(p),P30=float(tgres[k][-1,3]),bet30=float(tgres[k][-1,1:4].sum())) for k,p in tg.items()},
  tgm_thr_cut=1-1/R0_tgm, unc=unc)
json.dump(out,open(OUT/'summary.json','w'),indent=1)
print('wrote', OUT)
