#!/usr/bin/env python3
"""
TOX-102 (3-month B6C3F1 mouse, gavage) — Case 2 demo pipeline.
Covariate adjustment of organ (liver) weight by body weight.

INPUT (NTP CEBS individual-animal file):
  2032304_Individual_Animal_Organ_Weight_Data.xlsx   (sheet 'Data')
    columns: Concentration(dose mg/kg), Animal ID, Sex, Removal Day,
             Heart, Kidney-Right, Liver, Lungs, Terminal Body Weight, Testis Right, Thymus

OUTPUTS:
  tox102_3mo_mouse_organ_weights.csv
  tox102_case2_three_normalizations.png
  tox102_case2_covariate_scatter.png

Deps: pandas, openpyxl, numpy, scipy, statsmodels, matplotlib
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.stats import dunnett
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DOSES=[0,156,312,625,1250,2500]

# ---------- 1. PARSE + CLEAN ----------
def build_csv(path='2032304_Individual_Animal_Organ_Weight_Data.xlsx'):
    df=pd.read_excel(path,sheet_name='Data',header=0)
    df.columns=[str(c).strip() for c in df.columns]
    df=df.rename(columns={'Concentration':'Dose_mgkg','Animal ID':'Animal',
        'Terminal Body Weight':'TBW_g','Kidney-Right':'Kidney_R_g','Testis Right':'Testis_R_g',
        'Heart':'Heart_g','Liver':'Liver_g','Lungs':'Lungs_g','Thymus':'Thymus_g'})
    for c in ['Dose_mgkg','Heart_g','Kidney_R_g','Liver_g','Lungs_g','TBW_g','Testis_R_g','Thymus_g']:
        df[c]=pd.to_numeric(df[c],errors='coerce')
    df['Sex']=df['Sex'].map({'Male':'M','Female':'F'})
    df['RelLiver_pct']=(df['Liver_g']/df['TBW_g']*100).round(4)
    df['is_control']=df['Dose_mgkg'].eq(0)
    out=df[['Animal','Sex','Dose_mgkg','is_control','TBW_g','Liver_g','RelLiver_pct',
            'Heart_g','Kidney_R_g','Lungs_g','Thymus_g','Testis_R_g','Removal Day']]
    out=out.sort_values(['Sex','Dose_mgkg','Animal']).reset_index(drop=True)
    out.to_csv('tox102_3mo_mouse_organ_weights.csv',index=False)
    return out

# ---------- 2. STATS ----------
def dunnett_p(df,col,sex):
    sub=df[df.Sex==sex]; ctrl=sub[sub.Dose_mgkg==0][col].values
    res=dunnett(*[sub[sub.Dose_mgkg==d][col].values for d in DOSES[1:]],control=ctrl)
    return dict(zip(DOSES[1:],res.pvalue)), ctrl.mean()

def ancova(df,sex):
    sub=df[df.Sex==sex].copy(); sub['Dose']=sub['Dose_mgkg'].astype('category')
    full=smf.ols('Liver_g ~ C(Dose) + TBW_g',data=sub).fit()
    red =smf.ols('Liver_g ~ TBW_g',data=sub).fit()
    lr=anova_lm(red,full)
    return full.params['TBW_g'], full.pvalues['TBW_g'], lr['F'][1], lr['Pr(>F)'][1]

def report(df):
    for sex in ['M','F']:
        print('='*60); print({'M':'MALE','F':'FEMALE'}[sex]); print('='*60)
        for col,lab in [('Liver_g','ABSOLUTE liver (g)'),('RelLiver_pct','RELATIVE liver (%BW)')]:
            pmap,cm=dunnett_p(df,col,sex); sub=df[df.Sex==sex]
            print(f'\n{lab}: control={cm:.3f}')
            for d in DOSES[1:]:
                m=sub[sub.Dose_mgkg==d][col].mean(); p=pmap[d]
                s='***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
                print(f'  {d:>4}: {m:.3f} ({(m-cm)/cm*100:+.1f}%) p={p:.4f} {s}')
        slope,ps,F,pd_=ancova(df,sex)
        print(f'\nANCOVA Liver~Dose+BW: BW slope={slope:.4f} (p={ps:.2g}); '
              f'dose-after-BW F={F:.2f} p={pd_:.2g}')

# ---------- 3. FIGURES ----------
def figures(df):
    plt.rcParams.update({'font.size':11,'axes.edgecolor':'#444','axes.linewidth':.8,
        'savefig.dpi':160,'font.family':'DejaVu Sans'})
    C_ABS,C_REL,C_CTRL,GRID='#3a7d44','#b5651d','#9aa0a6','#e6e8eb'
    xpos=np.arange(len(DOSES)); jit=lambda n:(np.random.RandomState(3).rand(n)-.5)*.24

    # FIG 1: female, three normalizations
    fig,axes=plt.subplots(1,3,figsize=(15.5,5.2))
    sub=df[df.Sex=='F']
    def boxpanel(ax,col,color,title,pmap,cm,unit):
        for d,x in zip(DOSES,xpos):
            v=sub[sub.Dose_mgkg==d][col].values
            ax.boxplot(v,positions=[x],widths=.55,patch_artist=True,
                boxprops=dict(facecolor='#e8f0e9' if d==0 else 'white',color=color),
                medianprops=dict(color=color,lw=1.5),whiskerprops=dict(color=color),
                capprops=dict(color=color),showfliers=False,zorder=2)
            ax.scatter(x+jit(len(v)),v,s=20,color=color,alpha=.8,zorder=3,edgecolor='white',lw=.4)
            if d in pmap:
                p=pmap[d]; s='***' if p<.001 else '**' if p<.01 else '*' if p<.05 else ''
                if s: ax.text(x,v.max(),s,ha='center',va='bottom',fontsize=13,color=color,fontweight='bold')
        ax.axhline(cm,color=C_CTRL,ls='--',lw=.9,zorder=1)
        ax.set_title(title,fontsize=11,loc='left'); ax.set_xticks(xpos)
        ax.set_xticklabels([str(d) for d in DOSES],rotation=0); ax.set_xlabel('Dose (mg/kg/day)')
        ax.set_ylabel(unit); ax.yaxis.grid(True,color=GRID); ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
    pA,cA=dunnett_p(df,'Liver_g','F'); pR,cR=dunnett_p(df,'RelLiver_pct','F')
    boxpanel(axes[0],'Liver_g',C_ABS,'A. ABSOLUTE liver weight\n-> all doses n.s. ("no effect")',pA,cA,'Liver (g)')
    boxpanel(axes[1],'RelLiver_pct',C_REL,'B. RELATIVE liver (% body wt)\n-> significant from 625 ("clear effect")',pR,cR,'Liver / body weight (%)')
    # panel C: ANCOVA scatter
    ax=axes[2]
    cols=plt.cm.viridis(np.linspace(0,.9,len(DOSES)))
    for d,c in zip(DOSES,cols):
        s=sub[sub.Dose_mgkg==d]
        ax.scatter(s['TBW_g'],s['Liver_g'],s=30,color=c,label=str(d),edgecolor='white',lw=.4,zorder=3)
    # overall regression line
    z=np.polyfit(sub['TBW_g'],sub['Liver_g'],1); xs=np.linspace(sub['TBW_g'].min(),sub['TBW_g'].max(),50)
    ax.plot(xs,np.polyval(z,xs),color='#444',lw=1.3,ls='-',zorder=2)
    slope,ps,F,pd_=ancova(df,'F')
    ax.set_title(f'C. ANCOVA: liver vs body weight\n-> dose effect holds after BW adj. (p={pd_:.1e})',fontsize=11,loc='left')
    ax.set_xlabel('Terminal body weight (g)'); ax.set_ylabel('Liver (g)')
    ax.legend(title='Dose',fontsize=8,title_fontsize=8,frameon=False,loc='upper left',ncol=2)
    ax.yaxis.grid(True,color=GRID); ax.xaxis.grid(True,color=GRID); ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    fig.suptitle('TOX-102 - 3mo mouse (FEMALES) - same liver data, three normalizations, three conclusions',
                 fontsize=12.5,fontweight='bold',x=.01,ha='left')
    fig.tight_layout(rect=[0,0,1,.95])
    fig.savefig('tox102_case2_three_normalizations.png',bbox_inches='tight')

    # FIG 2: covariate scatter both sexes
    fig2,ax=plt.subplots(figsize=(8,5.4))
    for sex,mk,lab in [('M','o','Males'),('F','^','Females')]:
        s=df[df.Sex==sex]
        sc=ax.scatter(s['TBW_g'],s['Liver_g'],c=s['Dose_mgkg'],cmap='viridis',
                      marker=mk,s=34,edgecolor='white',lw=.4,label=lab)
    cb=plt.colorbar(sc,ax=ax); cb.set_label('Dose (mg/kg/day)')
    for sex,ls in [('M','-'),('F','--')]:
        s=df[df.Sex==sex]; z=np.polyfit(s['TBW_g'],s['Liver_g'],1)
        xs=np.linspace(s['TBW_g'].min(),s['TBW_g'].max(),50)
        ax.plot(xs,np.polyval(z,xs),color='#444',lw=1.2,ls=ls)
    ax.set_xlabel('Terminal body weight (g)'); ax.set_ylabel('Liver (g)')
    ax.set_title('Liver weight scales with body weight (slope ~0.043 g/g, p<1e-10):\nbody weight is a real covariate that must be adjusted for',fontsize=11,loc='left')
    ax.legend(frameon=False,loc='lower right')
    ax.yaxis.grid(True,color=GRID); ax.xaxis.grid(True,color=GRID); ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    fig2.tight_layout(); fig2.savefig('tox102_case2_covariate_scatter.png',bbox_inches='tight')

if __name__=='__main__':
    df=build_csv(); report(df); figures(df)
    print('\nDone: CSV + 2 figures written.')
