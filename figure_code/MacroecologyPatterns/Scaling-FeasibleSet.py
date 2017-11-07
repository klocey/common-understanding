from __future__ import division
import  matplotlib.pyplot as plt
from os.path import expanduser
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

def e_simpson(sad):
    sad = filter(lambda a: a != 0, sad)
    D = 0.0
    N = sum(sad)
    S = len(sad)
    for x in sad: D += (x*x) / (N*N)
    E = round((1.0/D)/S, 4)
    return E


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
        lab = r'$R_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$R_{microbes}$'+ ' = 2.34*'+r'$N$'+'$^{0.14}$'+'\n'
        lab += r'$R_{macrobes}$'+ ' = 1.7*'+r'$N$'+'$^{0.11}$'
        plt.text(0.2, 0.8, lab, fontsize=7)

    elif n == 2:
        lab = r'$D_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$D_{microbes}$'+ ' = 0.44*'+r'$N$'+'$^{0.92}$'+'\n'
        lab += r'$D_{macrobes}$'+ ' = 0.23*'+r'$N$'+'$^{0.99}$'
        plt.text(0.2, 3.0, lab, fontsize=7)

    elif n == 3:
        lab = r'$E_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$E_{microbes}$'+ ' = 0.58*'+r'$N$'+'$^{-0.23}$'+'\n'
        lab += r'$E_{macrobes}$'+ ' = 1.15*'+r'$N$'+'$^{-0.21}$'
        plt.text(0.2, -1.7, lab, fontsize=7)

    elif n == 4:
        lab = r'$S_{models}$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        lab += r'$S_{microbes}$'+ ' = 1.77*'+r'$N$'+'$^{0.38}$'+'\n'
        lab += r'$S_{macrobes}$'+ ' = 1.77*'+r'$N$'+'$^{0.24}$'
        plt.text(0.2, 1.9, lab, fontsize=7)

    #plt.hexbin(x2, y2, mincnt=1, gridsize = 40, bins='log', cmap=plt.cm.jet)
    plt.scatter(x2, y2, color = 'SkyBlue', alpha= 1 , s = 12, linewidths=0.5, edgecolor='Steelblue')

    if n == 3: plt.legend(loc='best', fontsize=6, frameon=False)

    plt.plot(x2, fitted,  color='k', ls='--', lw=1.0, alpha=0.9)
    plt.xlabel(xlab, fontsize=8)
    plt.ylabel(ylab, fontsize=8)
    plt.tick_params(axis='both', labelsize=5)
    plt.xlim(0, 1.05*max(x2))

    if n == 1: plt.ylim(0.0, 1.1)
    elif n == 2: plt.ylim(0.0, 4.2)
    elif n == 3: plt.ylim(-1.8, 0.05)
    elif n == 4: plt.ylim(0.4, 2.5)

    return fig


GenPath = mydir + "/results/simulated_data/"
dat = open(GenPath + 'partitions.csv', 'r')

Ns, Rs, Es, Ds, Ss = [], [], [], [], []
for line in dat:
    partition = eval(line)

    N = np.log10(sum(partition))
    D = np.log10(max(partition))
    S = np.log10(len(partition))
    E = np.log10(e_simpson(partition))
    R = stats.skew(partition)
    R = np.log10(abs(float(R)) + 1)
    if R < 0: R = R * -1

    Ns.append(N)
    Ds.append(D)
    Ss.append(S)
    Es.append(E)
    Rs.append(R)


fig = plt.figure()
xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Rarity ($R$), '+r'$log_{10}$'
fig = figplot(Ns, Rs, xlab, ylab, fig, 1)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Dominance ($D$), '+r'$log_{10}$'
fig = figplot(Ns, Ds, xlab, ylab, fig, 2)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Evenness ($E$), ' +r'$log_{10}$'
fig = figplot(Ns, Es, xlab, ylab, fig, 3)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = 'Richness ($S$), ' +r'$log_{10}$'
fig = figplot(Ns, Ss, xlab, ylab, fig, 4)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Partitions.png',
    dpi=400, bbox_inches = "tight")
plt.close()
