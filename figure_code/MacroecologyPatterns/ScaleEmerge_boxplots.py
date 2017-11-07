from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from math import floor
import scipy as sc
from scipy import stats
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def setBoxColors(bp, clr):
    plt.setp(bp['boxes'], color=clr)
    plt.setp(bp['caps'], color=clr)
    plt.setp(bp['whiskers'], color=clr)
    plt.setp(bp['medians'], color=clr)
    if clr == 'red': colors = ['red']*len(bp['boxes'])
    elif clr == 'blue': colors = ['blue']*len(bp['boxes'])
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.2)


m_points = [0.96, 0.45, 0.21, -0.35]
b_points = [0.29, 2.51, 1.15, 2.04]

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df2 = pd.DataFrame({'sim' : df['sim']})
df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).mean())
df2['D'] = np.log10(df['N.max'].groupby(df['sim']).mean())
df2['S'] = np.log10(df['species.richness'].groupby(df['sim']).mean())
df2['E'] = np.log10(df['simpson.e'].groupby(df['sim']).mean())
df2['R'] = df['logmod.skew'].groupby(df['sim']).mean()
df2 = df2[df2['S'] > 0]

def parse(x, y):
    mlist = []
    blist = []

    n1 = len(x)
    n2 = n1%100
    ns = [100]*int(floor(n1/100))
    ns.append(n2)
    print len(ns)
    for n in ns:
        x1 = x[:n]
        x = x[n:]
        y1 = y[:n]
        y = y[n:]
        m, b, r, p, std_err = stats.linregress(x1, y1)
        mlist.append(m)
        blist.append(b)
    return mlist, blist

R_slope, R_int = parse(df2['N'], df2['R'])
D_slope, D_int = parse(df2['N'], df2['D'])
E_slope, E_int = parse(df2['N'], df2['E'])
S_slope, S_int = parse(df2['N'], df2['S'])

data_to_plot1 = [D_slope, S_slope, R_slope, E_slope]
data_to_plot2 = [D_int,   S_int,   R_int,   E_int]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
bp = ax.boxplot(data_to_plot1, showfliers=False, patch_artist=True)
ax.set_xticklabels(['Dominance', 'Richness', 'Rarity', 'Evenness'], rotation=45)
setBoxColors(bp, 'red')
plt.ylabel('Scaling exponent')

for i, val in enumerate(m_points):
    plt.plot(i+1, val, 'k.', alpha=0.8, markersize=10)

'''
ax = fig.add_subplot(2, 1, 2)
bp = ax.boxplot(data_to_plot2, showfliers=False, patch_artist=True)
ax.set_xticklabels(['Rarity', 'Dominance', 'Evenness', 'Richness'], rotation=45)
setBoxColors(bp, 'blue')
plt.ylabel('Intercept')

for i, val in enumerate(b_points):
    plt.plot(i+1, val, 'k.', alpha=0.8, markersize=10)
'''

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/BoxPlots.png', dpi=600, bbox_inches = "tight")
plt.close()
