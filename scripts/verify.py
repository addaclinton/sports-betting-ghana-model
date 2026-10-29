"""Check the numbers in data/summary.json (and, optionally, those quoted in a manuscript .tex file).

Usage:  python scripts/verify.py                      # numerical self-checks only
        python scripts/verify.py manuscript.tex       # also check numbers/citations in the manuscript
"""
import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TEX=sys.argv[1] if len(sys.argv)>1 else None
s=open(TEX).read() if TEX else ""
J=json.load(open(ROOT/'data'/'summary.json'))
errs=[]
def need(label,cond): 
    if not cond: errs.append(label)
    print(("OK   " if cond else "FAIL ")+label)
def has(x): return (x in s) if TEX else True
# ---- citations vs bibliography
cites=set(k for m in re.findall(r'\\cite\{([^}]*)\}',s) for k in m.split(','))
bibs=set(re.findall(r'\\bibitem\{([^}]*)\}',s))
need("every cited key has a bibitem: missing=%s"%sorted(cites-bibs), (not (cites-bibs)) or not TEX)
need("no uncited bibitem: unused=%s"%sorted(bibs-cites), (not (bibs-cites)) or not TEX)
# ---- numbers
u=J['unc']; sc=J['scen']; tg=J['tgm_scen']; e=J['eq']; en=J['end30']; N0=17931673
f0=lambda x:f"{x:,.0f}"
need("R0 base 1.30", abs(J['R0']-1.3)<1e-9)
need("beta 0.2405",  has("0.2405") and round(J['beta'],4)==0.2405)
need("alpha1 0.6351, alpha2 0.1100", round(J['a1'],4)==0.6351 and round(J['a2'],4)==0.1100)
need("mu 0.02097, Lambda 376,005", round(J['mu'],5)==0.02097 and f0(J['Lam'])=="376,005")
need("mu = 1/47.69", abs(1/47.69-J['mu'])<2e-5)
for nm,v in zip(['S*','E*','R*','P*','Q*'],e): print("   eq",nm,f0(v))
need("E* ~337,900", f0(round(e[1],-2))=="337,900" and has("337{,}900"))
need("R* ~929,300", f0(round(e[2],-2))=="929,300" and has("929{,}300"))
need("P* ~844,800", f0(round(e[3],-2))=="844,800" and has("844{,}800"))
need("Q* ~2.03 million", round(e[4]/1e6,2)==2.03 and has("2.03$ million"))
need("bettor prev 11.8%", round(100*sum(e[1:4])/N0,1)==11.8 and has("11.8\\%"))
need("problem share 4.7%", round(100*e[3]/N0,1)==4.7 and has("4.7\\%"))
need("pop constant", abs(sum(e)-N0)<1)
need("30y E 104,800",  round(en[1],-2)==104800 and has("104,800"))
need("30y R 256,300",  round(en[2],-2)==256300 and has("256,300"))
need("30y P 186,000",  round(en[3],-3)==186000 and has("186,000"))
need("30y Q 480,600",  round(en[4],-2)==480600 and has("480,600"))
need("30y bettors 547,100 / 3.1%", round(sum(en[1:4]),-2)==547100 and round(100*sum(en[1:4])/N0,1)==3.1 and has("547,100"))
need("a=1.16", round(J['a1']+0.5+J['mu'],2)==1.16)
sn=J['sens']
for k,v in [('beta','+1.0000'),('a1','+0.2906'),('a2','0.0000'),('g1','$-0.4325$'),('g2','$-0.3637$'),('g3','$-0.3307$'),('mu','$-0.1637$')]:
    vv=sn[k]; txt=(f"{vv:+.4f}" if vv>1e-6 else (f"$-{abs(vv):.4f}$" if vv<-1e-6 else "0.0000"))
    need(f"sensitivity {k}={txt}", has(v) and (txt==v))
# interventions
red=lambda a,b:100*(1-a/b)
need("beta30 P30 57,224, R0 .910, red 69.2", f0(sc['beta30']['P30'])=="57,224" and round(sc['beta30']['R0'],3)==0.910 and round(red(sc['beta30']['P30'],sc['base']['P30']),1)==69.2)
need("g3 P30 109,065, R0 1.148, red 41.4", f0(sc['g3_50']['P30'])=="109,065" and round(sc['g3_50']['R0'],3)==1.148 and round(red(sc['g3_50']['P30'],sc['base']['P30']),1)==41.4)
need("base P30 186,025", f0(sc['base']['P30'])=="186,025")
need("peak ~68,000 at yr 13", round(J['b30_peak'],-3)==68000 and round(J['b30_peak_t'])==13)
need("beta cut needed 23.1%", round(100*(1-1/1.3),1)==23.1 and has("23.1\\%"))
need("g3 for R0=1: 0.265, 2.65-fold", round(J['g3_for_R0_1'],3)==0.265 and round(J['g3_for_R0_1']/0.1,2)==2.65)
need("TGM R0 5.49, beta 1.015", round(J['R0_tgm'],2)==5.49 and round(J['tgm_beta'],3)==1.015)
need("TGM shares: P 16.7%, Q 40.0%", round(100*J['tgm_eq_share'][3],1)==16.7 and round(100*J['tgm_eq_share'][4],1)==40.0)
need("TGM bettors = 41.74%", round(100*sum(J['tgm_eq_share'][1:4]),2)==41.74)
need("TGM scen base 3,008,262", f0(tg['base']['P30'])=="3,008,262")
need("TGM beta30 2,728,026 R0 3.842 red 9.3", f0(tg['beta30']['P30'])=="2,728,026" and round(tg['beta30']['R0'],3)==3.842 and round(red(tg['beta30']['P30'],tg['base']['P30']),1)==9.3)
need("TGM g3 2,188,505 R0 4.846 red 27.3", f0(tg['g3_50']['P30'])=="2,188,505" and round(tg['g3_50']['R0'],3)==4.846 and round(red(tg['g3_50']['P30'],tg['base']['P30']),1)==27.3)
need("TGM needed cut 81.8%", round(100*J['tgm_thr_cut'],1)==81.8 and has("81.8\\%"))
need("share R0>1 = 88.9%", round(100*u['share_gt1'],1)==88.9)
need("P30 median 1.23m, 5-95: 0.09-3.74m", round(u['P30_med']/1e6,2)==1.23 and round(u['P30_lo']/1e6,2)==0.09 and round(u['P30_hi']/1e6,2)==3.74)
need("median 6.8% of adults", round(100*u['P30_med']/N0,1)==6.8)
pr=u['prcc']
need("PRCC .92 .54 .36 -.33 .00", [round(pr['R0'],2),round(pr['delta'],2),round(pr['g3'],2),round(pr['g1'],2),round(pr['g2'],2)]==[0.92,0.54,0.36,-0.33,0.0])
need("beta reduction 52% (9-88)", [round(100*u['redb_med']),round(100*u['redb_lo']),round(100*u['redb_hi'])]==[52,9,88])
need("g3 reduction 31% (13-62)", [round(100*u['redg_med']),round(100*u['redg_lo']),round(100*u['redg_hi'])]==[31,13,62])
need("68% beta better", round(100*u['frac_beta_better'])==68)
need("R0<1 in 19.9% vs 6.7%", round(100*u['below1_beta'],1)==19.9 and round(100*u['below1_g3'],1)==6.7)
need("vector vs scalar < 1e-10", u['check_vs_scalar']<1e-10)
need("log conversions", round(-__import__('math').log(1-0.416),2)==0.54 and round(-__import__('math').log(0.76)/8,3)==0.034 and round(-__import__('math').log(1-0.584),2)==0.88)
need("58.4% = 29.2+29.2; 41.6+58.4=100", abs(29.2+29.2-58.4)<1e-9 and abs(41.6+58.4-100)<1e-9)
need("pop share 58.2%: 17,931,673/30,832,019", round(100*17931673/30832019,1)==58.2)
need("composition 16:44:40 sums 100; 2.75 and 0.40/0.44", 16+44+40==100 and 44/16==2.75)
print("\nERRORS:",errs)
sys.exit(1 if errs else 0)
