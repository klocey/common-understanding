from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
import scipy as sc
from scipy import stats
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def figplot(x, y, sims, clrs, xlab, ylab, fig, n):

    fig.add_subplot(2, 2, n)
    y2 = list(y)
    x2 = list(x)
    clrs = list(clrs)

    plt.scatter(x2, y2, color = clrs, s = 2, linewidths=0.0)

    d = pd.DataFrame({'x': list(x2)})
    d['y'] = list(y2)
    f = smf.ols('y ~ x', d).fit()
    st, data, ss2 = summary_table(f, alpha=0.05)
    fitted = data[:,2]
    m, b, r, p, std_err = stats.linregress(x2, y2)

    if n == 1:
        lab = r'$R_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$R_{microbes}$'+ ' = 2.34*'+r'$N$'+'$^{0.14}$'+'\n'
        lab += r'$R_{macrobes}$'+ ' = 1.7*'+r'$N$'+'$^{0.11}$'
        plt.text(0.2, 1.4, lab, fontsize=7)

    elif n == 2:
        lab = r'$D_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$D_{microbes}$'+ ' = 0.44*'+r'$N$'+'$^{0.92}$'+'\n'
        lab += r'$D_{macrobes}$'+ ' = 0.23*'+r'$N$'+'$^{0.99}$'
        plt.text(0.2, 2.5, lab, fontsize=7)

    elif n == 3:
        lab = r'$E_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$E_{microbes}$'+ ' = 0.58*'+r'$N$'+'$^{-0.23}$'+'\n'
        lab += r'$E_{macrobes}$'+ ' = 1.15*'+r'$N$'+'$^{-0.21}$'
        plt.text(0.2, -3.4, lab, fontsize=7)

    elif n == 4:
        lab = r'$S_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$S_{microbes}$'+ ' = 1.77*'+r'$N$'+'$^{0.38}$'+'\n'
        lab += r'$S_{macrobes}$'+ ' = 1.77*'+r'$N$'+'$^{0.24}$'
        plt.text(0.2, 2.8, lab, fontsize=7)

    if n == 3: plt.legend(loc='best', fontsize=6, frameon=False)

    plt.plot(x2, fitted,  color='k', ls='--', lw=1.0, alpha=0.9)
    plt.xlabel(xlab, fontsize=8)
    plt.ylabel(ylab, fontsize=8)
    plt.tick_params(axis='both', labelsize=5)
    plt.xlim(0, 1.05*max(x2))

    if n == 1: plt.ylim(0.0, max(y2))
    elif n == 2: plt.ylim(0.0, max(y2))
    elif n == 3: plt.ylim(min(y2), 0)
    elif n == 4: plt.ylim(0.4, max(y2))

    return fig

df = pd.read_csv(mydir + '/results/simulated_data/SimData-Random.csv')
df = df[df['species.richness'] > 0]
df = df[df['total.abundance'] > 0]
df = df[df['logmod.skew'] != 0]
df = df[df['simpson.e'] != 0]
#df = df[df['sim'] < 10]

df2 = pd.DataFrame({'sim' : df['sim']})
df2['N'] = np.log10(df['total.abundance'])
df2['D'] = np.log10(df['N.max'])
df2['S'] = np.log10(df['species.richness'])
df2['E'] = np.log10(df['simpson.e'])
df2['R'] = df['logmod.skew']
df2['sim'] = df['sim']
df2['frac'] = df['frac']
df2['choose'] = df['choose']
df2['color'] = [s.replace('\'', '') for s in df['clr']]

df3 = df2.replace([np.inf, -np.inf], np.nan).dropna()

fig = plt.figure()
xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Rarity ($R$), '+r'$log_{10}$'
fig = figplot(df3['N'], df3['R'], df3['sim'], df3['color'], xlab, ylab, fig, 1)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Dominance ($D$), '+r'$log_{10}$'
fig = figplot(df3['N'], df3['D'], df3['sim'], df3['color'], xlab, ylab, fig, 2)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Evenness ($E$), ' +r'$log_{10}$'
fig = figplot(df3['N'], df3['E'], df3['sim'], df3['color'], xlab, ylab, fig, 3)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Richness ($S$), ' +r'$log_{10}$'
fig = figplot(df3['N'], df3['S'], df3['sim'], df3['color'], xlab, ylab, fig, 4)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Random.png',
    dpi=400, bbox_inches = "tight")
plt.close()
