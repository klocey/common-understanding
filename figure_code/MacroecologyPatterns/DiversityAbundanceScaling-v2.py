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


def xfrm(X, _max): return -np.log10(_max - np.array(X))

def figplot(x, y, xlab, ylab, fig, n):

    '''main figure plotting function'''

    fig.add_subplot(2, 2, n)
    y2 = list(y)
    x2 = list(x)

    d = pd.DataFrame({'x': list(x2)})
    d['y'] = list(y2)
    f = smf.ols('y ~ x', d).fit()

    m, b, r, p, std_err = stats.linregress(x2, y2)

    st, data, ss2 = summary_table(f, alpha=0.05)
    fitted = data[:,2]
    mean_ci_low, mean_ci_upp = data[:,4:6].T
    ci_low, ci_upp = data[:,6:8].T

    x2, y2, fitted, ci_low, ci_upp = zip(*sorted(zip(x2, y2, fitted, ci_low, ci_upp)))

    if n == 1:
        lab = r'$rarity$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$micro$'+ ' = 2.34*'+r'$N$'+'$^{0.14}$'+'\n'
        lab += r'$macro$'+ ' = 1.7*'+r'$N$'+'$^{0.11}$'
        plt.text(0.2, 0.8, lab, fontsize=8)

    elif n == 2:
        lab = r'$Nmax$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$micro$'+ ' = 0.44*'+r'$N$'+'$^{0.92}$'+'\n'
        lab += r'$macro$'+ ' = 0.23*'+r'$N$'+'$^{0.99}$'
        plt.text(0.2, 3.0, lab, fontsize=8)

    elif n == 3:
        lab = r'$Ev$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$micro$'+ ' = 0.58*'+r'$N$'+'$^{-0.23}$'+'\n'
        lab += r'$macro$'+ ' = 1.15*'+r'$N$'+'$^{-0.21}$'
        plt.text(0.2, -1.55, lab, fontsize=8)

    elif n == 4:
        lab = r'$S$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$micro$'+ ' = 1.77*'+r'$N$'+'$^{0.38}$'+'\n'
        lab += r'$macro$'+ ' = 1.77*'+r'$N$'+'$^{0.24}$'
        plt.text(0.2, 1.8, lab, fontsize=8)

    #plt.scatter(x2, y2, color = 'SkyBlue', alpha= 1 , s = 12, linewidths=0.5, edgecolor='Steelblue')
    plt.hexbin(x2, y2, mincnt=1, gridsize = 40, bins='log', cmap=plt.cm.jet)

    if n == 3: plt.legend(loc='best', fontsize=6, frameon=False)

    #plt.fill_between(x2, ci_upp, ci_low, color='b', lw=0.1, alpha=0.15)
    plt.plot(x2, fitted,  color='b', ls='--', lw=1.0, alpha=0.9)
    plt.xlabel(xlab, fontsize=8)
    plt.ylabel(ylab, fontsize=8)
    plt.tick_params(axis='both', labelsize=5)
    plt.xlim(0, 1.05*max(x2))

    if n == 1: plt.ylim(0.0, 1.1)
    elif n == 2: plt.ylim(0.0, 4.2)
    elif n == 3: plt.ylim(-1.8, 0.05)
    elif n == 4: plt.ylim(0.4, 2.5)

    return fig

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df2 = pd.DataFrame({'sim' : df['sim']})

df = df[df['species.richness'] > 0]
df = df[df['total.abundance'] > 0]
df = df[df['logmod.skew'] != 0]
df = df[df['simpson.e'] != 0]

df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).mean())
df2['D'] = np.log10(df['N.max'].groupby(df['sim']).mean())
df2['S'] = np.log10(df['species.richness'].groupby(df['sim']).mean())
df2['E'] = np.log10(df['simpson.e'].groupby(df['sim']).mean())
df2['R'] = df['logmod.skew'].groupby(df['sim']).mean()
df2['gr'] = df['gr'].groupby(df['sim']).mean()
df2['mt'] = df['mt'].groupby(df['sim']).mean()
df2['q'] = df['q'].groupby(df['sim']).mean()

df2 = df2[df2['q'] != 3]
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

fig = plt.figure()
xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Rarity, '+r'$log_{10}$'
fig = figplot(df2['N'], df2['R'], xlab, ylab, fig, 1)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Dominance, '+r'$log_{10}$'
fig = figplot(df2['N'], df2['D'], xlab, ylab, fig, 2)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Evenness, ' +r'$log_{10}$'
fig = figplot(df2['N'], df2['E'], xlab, ylab, fig, 3)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Richness, ' +r'$log_{10}$'
fig = figplot(df2['N'], df2['S'], xlab, ylab, fig, 4)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling.png',
    dpi=400, bbox_inches = "tight")
plt.close()
