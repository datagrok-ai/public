#!/usr/bin/env python3
"""
TR-594 (3-month rat, inhalation) — Case 1 demo pipeline.
Builds a cleaned long-format terminal-body-weight CSV and the demo figures
illustrating why the control must be stratified by sex.

INPUTS (NTP CEBS individual-animal .xls files = HTML tables):
  1047201_Male_Individual_Animal_Body_Weight_Data.xls
  1047201_Female_Individual_Animal_Body_Weight_Data.xls

OUTPUTS:
  tr594_3mo_rat_terminal_bw.csv
  tr594_case1_naive_vs_correct.png
  tr594_case1_trends.png

Deps: pandas, lxml, numpy, scipy, matplotlib
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import dunnett, norm
import warnings; warnings.filterwarnings('ignore')

DOSE_MAP = {'Vehicle Control':0,'125ppm':125,'250ppm':250,
            '500ppm':500,'1000ppm':1000,'2000ppm':2000}
DOSES = [0,125,250,500,1000,2000]

# ---------- 1. PARSE + CLEAN ----------
def load_bw(path, sex):
    t = pd.read_html(path)[1]                 # table[1] = data; table[0] = metadata
    t.columns = [str(c).strip() for c in t.columns]
    return pd.DataFrame({
        'Animal':        t['Animal Number'].astype(int),
        'Sex':           sex,
        'Dose_ppm':      t['Dose Group'].map(DOSE_MAP).astype(int),
        'Dose_label':    t['Dose Group'],
        'Days_on_study': t['Days On Study'].astype(int),
        'Terminal_BW_g': pd.to_numeric(t['Week 14'], errors='coerce'),
        'Removal_reason':t['Removal Reason'],
    })

def build_csv():
    m = load_bw('1047201_Male_Individual_Animal_Body_Weight_Data.xls','M')
    f = load_bw('1047201_Female_Individual_Animal_Body_Weight_Data.xls','F')
    df = pd.concat([m,f], ignore_index=True)
    df['is_control'] = df['Dose_ppm'].eq(0)
    df = df.sort_values(['Sex','Dose_ppm','Animal']).reset_index(drop=True)
    df.to_csv('tr594_3mo_rat_terminal_bw.csv', index=False)
    return df

# ---------- 2. STATS ----------
def jonckheere(groups):
    """Jonckheere-Terpstra trend test. groups ordered by dose -> (z, two-sided p)."""
    k=len(groups); U=0
    for i in range(k):
        for j in range(i+1,k):
            for a in groups[i]:
                for b in groups[j]:
                    U += (b>a) + 0.5*(b==a)
    ns=[len(g) for g in groups]; N=sum(ns)
    EU=(N**2-sum(n**2 for n in ns))/4
    VU=(N**2*(2*N+3)-sum(n**2*(2*n+3) for n in ns))/72
    z=(U-EU)/np.sqrt(VU)
    return z, 2*(1-norm.cdf(abs(z)))

def analyze(ana, subset_label):
    ctrl = ana[ana.Dose_ppm==0]['Terminal_BW_g'].values
    treat=[ana[ana.Dose_ppm==d]['Terminal_BW_g'].values for d in DOSES[1:]]
    res=dunnett(*treat, control=ctrl)
    z,p=jonckheere([ana[ana.Dose_ppm==d]['Terminal_BW_g'].values for d in DOSES])
    print(f"\n[{subset_label}] control mean={ctrl.mean():.1f} SD={ctrl.std():.1f} N={len(ctrl)}")
    for d,pv in zip(DOSES[1:],res.pvalue):
        mn=ana[ana.Dose_ppm==d]['Terminal_BW_g'].mean()
        print(f"   {d:>4} ppm  mean={mn:5.1f}  d={mn-ctrl.mean():+6.1f}  Dunnett p={pv:.4f}")
    print(f"   Jonckheere z={z:+.2f} p={p:.4f}")
    return res

# ---------- 3. FIGURES ----------
def figures(df):
    plt.rcParams.update({'font.size':11,'axes.edgecolor':'#444','axes.linewidth':0.8,
        'savefig.dpi':160,'font.family':'DejaVu Sans'})
    C_M,C_F,C_CTRL,GRID='#2a7fb8','#d4682a','#9aa0a6','#e6e8eb'
    xpos=np.arange(len(DOSES))
    jit=lambda n: (np.random.RandomState(7).rand(n)-0.5)*0.24

    fig,axes=plt.subplots(1,3,figsize=(15,5.2))
    ax=axes[0]
    for d,x in zip(DOSES,xpos):
        for sex,col in [('M',C_M),('F',C_F)]:
            v=df[(df.Dose_ppm==d)&(df.Sex==sex)]['Terminal_BW_g'].values
            ax.scatter(x+jit(len(v)),v,s=22,color=col,alpha=.75,zorder=3,edgecolor='white',linewidth=.4)
        allv=df[df.Dose_ppm==d]['Terminal_BW_g'].values
        ax.boxplot(allv,positions=[x],widths=.5,patch_artist=True,
            boxprops=dict(facecolor='none',color='#555'),medianprops=dict(color='#222'),
            whiskerprops=dict(color='#888'),capprops=dict(color='#888'),showfliers=False,zorder=2)
    ax.set_title('A. Naive: sexes pooled\n(control bimodal -> effect cancels)',fontsize=11,loc='left')
    ax.text(.02,.04,'Jonckheere p = 0.33 (no trend)\ncontrol SD = 79 g',transform=ax.transAxes,
            fontsize=9.5,color='#b00',va='bottom')

    def panel(ax,sex,col,title,jp,arrow):
        sub=df[df.Sex==sex]
        ctrl=sub[sub.Dose_ppm==0]['Terminal_BW_g'].values
        res=dunnett(*[sub[sub.Dose_ppm==d]['Terminal_BW_g'].values for d in DOSES[1:]],control=ctrl)
        pmap=dict(zip(DOSES[1:],res.pvalue))
        for d,x in zip(DOSES,xpos):
            v=sub[sub.Dose_ppm==d]['Terminal_BW_g'].values
            ax.boxplot(v,positions=[x],widths=.55,patch_artist=True,
                boxprops=dict(facecolor='#dfe9f0' if d==0 else 'white',color=col),
                medianprops=dict(color=col,linewidth=1.5),whiskerprops=dict(color=col),
                capprops=dict(color=col),showfliers=False,zorder=2)
            ax.scatter(x+jit(len(v)),v,s=22,color=col,alpha=.8,zorder=3,edgecolor='white',linewidth=.4)
            if d in pmap:
                p=pmap[d]; s='***' if p<.001 else '**' if p<.01 else '*' if p<.05 else ''
                if s: ax.text(x,v.max()+6,s,ha='center',fontsize=13,color=col,fontweight='bold')
        ax.axhline(ctrl.mean(),color=C_CTRL,ls='--',lw=.9,zorder=1)
        ax.set_title(title,fontsize=11,loc='left')
        ax.text(.02,.96,f'Jonckheere p {jp}  {arrow}',transform=ax.transAxes,
                fontsize=9.5,color=col,va='top',fontweight='bold')
    panel(axes[1],'M',C_M,'B. Correct - males (matched control)','= 0.029','down')
    panel(axes[2],'F',C_F,'C. Correct - females (matched control)','< 0.0001','up')

    for ax in axes:
        ax.set_xticks(xpos); ax.set_xticklabels([str(d) for d in DOSES])
        ax.set_xlabel('Dose (ppm)'); ax.set_axisbelow(True)
        ax.yaxis.grid(True,color=GRID); ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Terminal body weight (g)')
    fig.legend([Line2D([],[],marker='o',ls='',color=C_M),Line2D([],[],marker='o',ls='',color=C_F)],
               ['Males','Females'],loc='upper right',frameon=False,ncol=2,bbox_to_anchor=(.99,1.0))
    fig.suptitle('TR-594 - 3mo rat - terminal body weight: why the control must be sex-stratified',
                 fontsize=12.5,fontweight='bold',x=.01,ha='left')
    fig.tight_layout(rect=[0,0,1,.95])
    fig.savefig('tr594_case1_naive_vs_correct.png',bbox_inches='tight')

    fig2,ax=plt.subplots(figsize=(7.5,5))
    for sex,col,lab in [('M',C_M,'Males'),('F',C_F,'Females')]:
        sub=df[df.Sex==sex]
        means=[sub[sub.Dose_ppm==d]['Terminal_BW_g'].mean() for d in DOSES]
        ses=[sub[sub.Dose_ppm==d]['Terminal_BW_g'].sem() for d in DOSES]
        ax.errorbar(xpos,means,yerr=ses,marker='o',color=col,capsize=3,lw=2,label=lab,zorder=3)
    pooled=[df[df.Dose_ppm==d]['Terminal_BW_g'].mean() for d in DOSES]
    ax.plot(xpos,pooled,color='#888',ls=':',lw=2,marker='s',label='Pooled (M+F)',zorder=2)
    ax.set_xticks(xpos); ax.set_xticklabels([str(d) for d in DOSES])
    ax.set_xlabel('Dose (ppm)'); ax.set_ylabel('Mean terminal BW (g) +/- SE')
    ax.set_title('Sex-divergent dose response; pooled line (dotted) is nearly flat',fontsize=11,loc='left')
    ax.legend(frameon=False); ax.yaxis.grid(True,color=GRID); ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    fig2.tight_layout(); fig2.savefig('tr594_case1_trends.png',bbox_inches='tight')

if __name__=='__main__':
    df=build_csv()
    ana=df.dropna(subset=['Terminal_BW_g'])         # drop moribund (no terminal weight)
    analyze(ana[ana.Sex=='M'],'MALES (correct)')
    analyze(ana[ana.Sex=='F'],'FEMALES (correct)')
    analyze(ana,'POOLED (naive)')
    figures(df)
    print('\nDone: CSV + 2 figures written.')
