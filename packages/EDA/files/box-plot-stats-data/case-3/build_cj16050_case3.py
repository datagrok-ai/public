#!/usr/bin/env python3
"""
CJ16050 (PHUSE SEND example, respiratory function in rats) — Case 3 demo pipeline.
Small-N display (strip vs box) + dose-response trend testing (Jonckheere-Terpstra)
vs pairwise control comparisons (Dunnett).

INPUTS (SEND .xpt domains, downloaded from PHUSE phuse-scripts GitHub):
  dm.xpt (demographics: SEX, SETCD), re.xpt (respiratory endpoints: RESTRESN, RETESTCD, RETPT)
  (ts/tx/ta/te/ex/se/cl/ds also exist in the study but are not needed here)

OUTPUTS:
  cj16050_resp_function.csv
  cj16050_case3_box_vs_strip.png
  cj16050_case3_trend_vs_pairwise.png

Deps: pandas, numpy, scipy, matplotlib  (pandas.read_sas reads XPORT natively)

Study structure: 3 dose groups (0 / 100 / 1000 mg/kg Compound A) x 6 MALES = 18 rats.
Endpoints: Respiratory Rate (breaths/min), Tidal Volume (mL/breath), Minute Volume (mL/min)
at Predose / 1 / 2 / 4 / 8 h postdose on Day 1. (No females -> sex faceting not applicable;
that knob is covered by Cases 1-2.)
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.stats import dunnett, norm
import matplotlib.pyplot as plt

DOSES=[0,100,1000]
DOSE_MAP={'00':0,'01':100,'02':1000}

# ---------- 1. PARSE SEND -> CSV ----------
def build_csv():
    dm=pd.read_sas('dm.xpt',format='xport',encoding='latin1')[['USUBJID','SEX','SETCD','ARM']]
    re=pd.read_sas('re.xpt',format='xport',encoding='latin1')
    dm['Dose_mgkg']=dm['SETCD'].astype(str).str.strip().map(DOSE_MAP)
    dm['SEX']=dm['SEX'].astype(str).str.strip()
    re=re.merge(dm,on='USUBJID')
    re['RESTRESN']=pd.to_numeric(re['RESTRESN'],errors='coerce')
    for c in ['RETPT','RETESTCD','RETEST','RESTRESU']: re[c]=re[c].astype(str).str.strip()
    out=pd.DataFrame({
        'USUBJID':re['USUBJID'].astype(str).str.strip(),
        'Animal':re['USUBJID'].astype(str).str.strip().str[-3:],
        'Sex':re['SEX'],'Dose_mgkg':re['Dose_mgkg'],
        'Dose_label':re['ARM'].astype(str).str.strip(),
        'TestCD':re['RETESTCD'],'Test':re['RETEST'],
        'Timepoint':re['RETPT'],'TimepointNum':pd.to_numeric(re['RETPTNUM'],errors='coerce'),
        'Value':re['RESTRESN'],'Unit':re['RESTRESU'].replace('nan',''),
        'is_control':re['Dose_mgkg'].eq(0)})
    out=out.dropna(subset=['Value']).sort_values(['TestCD','Timepoint','Dose_mgkg','Animal']).reset_index(drop=True)
    out.to_csv('cj16050_resp_function.csv',index=False)
    return out

# ---------- 2. STATS ----------
def jonckheere(groups):
    k=len(groups);U=0
    for i in range(k):
        for j in range(i+1,k):
            for a in groups[i]:
                for b in groups[j]:
                    U+=(b>a)+0.5*(b==a)
    ns=[len(g) for g in groups];N=sum(ns)
    EU=(N**2-sum(n*n for n in ns))/4
    VU=(N**2*(2*N+3)-sum(n*n*(2*n+3) for n in ns))/72
    z=(U-EU)/np.sqrt(VU);return z,2*(1-norm.cdf(abs(z)))

def grp(df,test,tp):
    sub=df[(df.TestCD==test)&(df.Timepoint==tp)]
    return [sub[sub.Dose_mgkg==d]['Value'].values for d in DOSES]

def stats(df,test,tp):
    g=grp(df,test,tp)
    res=dunnett(g[1],g[2],control=g[0])
    z,p=jonckheere(g)
    return g,dict(zip(DOSES[1:],res.pvalue)),(z,p)

# ---------- 3. FIGURES ----------
def figures(df):
    plt.rcParams.update({'font.size':11,'axes.edgecolor':'#444','axes.linewidth':.8,
        'savefig.dpi':160,'font.family':'DejaVu Sans'})
    GRID='#e6e8eb'; C1='#3a6ea5'; xs=np.arange(len(DOSES))
    jit=lambda n:(np.random.RandomState(5).rand(n)-.5)*.22
    TP='4-hour postdose'

    def stars(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'

    # FIG 1: box vs strip, Tidal Volume
    g,pm,(z,jp)=stats(df,'TIDALVOL',TP)
    fig,axes=plt.subplots(1,2,figsize=(11,5),sharey=True)
    # box
    ax=axes[0]
    ax.boxplot(g,positions=xs,widths=.5,patch_artist=True,
        boxprops=dict(facecolor='#eef2f6',color=C1),medianprops=dict(color=C1,lw=1.6),
        whiskerprops=dict(color=C1),capprops=dict(color=C1),showfliers=False)
    ax.set_title('A. Box plot, N=6/group\n(quartiles on 6 points imply false precision)',fontsize=11,loc='left')
    # strip
    ax=axes[1]
    for x,d,gg in zip(xs,DOSES,g):
        ax.scatter(x+jit(len(gg)),gg,s=42,color=C1,alpha=.85,edgecolor='white',lw=.5,zorder=3)
        ax.hlines(gg.mean(),x-.2,x+.2,color=C1,lw=2,zorder=2)
        if d in pm:
            s=stars(pm[d])
            ax.text(x,max(gg)+.03,s,ha='center',fontsize=12,color=C1,fontweight='bold')
    ax.set_title('B. Strip plot, same data\n(shows the real spread + every animal)',fontsize=11,loc='left')
    ax.text(.98,.96,f'Jonckheere p={jp:.3f}\nDunnett 1000 vs 0: {stars(pm[1000])}',
            transform=ax.transAxes,ha='right',va='top',fontsize=9.5,color=C1)
    for ax in axes:
        ax.set_xticks(xs); ax.set_xticklabels([str(d) for d in DOSES]); ax.set_xlabel('Dose (mg/kg)')
        ax.yaxis.grid(True,color=GRID); ax.set_axisbelow(True); ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Tidal volume (mL/breath)')
    fig.suptitle('CJ16050 - Tidal Volume @ 4h - small-N honesty: box vs strip',
                 fontsize=12.5,fontweight='bold',x=.01,ha='left')
    fig.tight_layout(rect=[0,0,1,.95]); fig.savefig('cj16050_case3_box_vs_strip.png',bbox_inches='tight')

    # FIG 2: monotonic (Jonckheere works) vs biphasic (Jonckheere fails, Dunnett works)
    fig2,axes=plt.subplots(1,2,figsize=(12,5))
    panels=[('TIDALVOL','Tidal volume (mL/breath)','A. MONOTONIC decrease\nJonckheere DETECTS the trend',C1),
            ('RESPRATE','Respiratory rate (breaths/min)','B. BIPHASIC (up then down)\nJonckheere MISSES it; Dunnett catches both','#b5651d')]
    for ax,(test,ylab,title,col) in zip(axes,panels):
        g,pm,(z,jp)=stats(df,test,TP)
        for x,d,gg in zip(xs,DOSES,g):
            ax.scatter(x+jit(len(gg)),gg,s=42,color=col,alpha=.85,edgecolor='white',lw=.5,zorder=3)
            ax.hlines(gg.mean(),x-.2,x+.2,color=col,lw=2,zorder=2)
            if d in pm:
                ax.text(x,max(gg)+(max(gg)-min(gg))*.05+.01,stars(pm[d]),ha='center',fontsize=12,color=col,fontweight='bold')
        ax.plot(xs,[gg.mean() for gg in g],color=col,lw=1.2,ls=':',zorder=1)
        ax.set_title(title,fontsize=11,loc='left'); ax.set_xticks(xs)
        ax.set_xticklabels([str(d) for d in DOSES]); ax.set_xlabel('Dose (mg/kg)'); ax.set_ylabel(ylab)
        ax.text(.02,.97,f'Jonckheere p={jp:.3f} ({stars(jp)})',transform=ax.transAxes,va='top',
                fontsize=9.5,color=col,fontweight='bold')
        ax.yaxis.grid(True,color=GRID); ax.set_axisbelow(True); ax.spines[['top','right']].set_visible(False)
    fig2.suptitle('CJ16050 @ 4h - a trend test is not a silver bullet: it only sees monotonic effects',
                  fontsize=12.5,fontweight='bold',x=.01,ha='left')
    fig2.tight_layout(rect=[0,0,1,.95]); fig2.savefig('cj16050_case3_trend_vs_pairwise.png',bbox_inches='tight')

if __name__=='__main__':
    df=build_csv()
    for test in ['TIDALVOL','RESPRATE','MV']:
        g,pm,(z,jp)=stats(df,test,'4-hour postdose')
        print(f"{test} @4h: means={[round(x.mean(),2) for x in g]} "
              f"Dunnett(100,1000)=({pm[100]:.3f},{pm[1000]:.3f}) Jonckheere p={jp:.4f}")
    figures(df)
    print('Done: CSV + 2 figures written.')
