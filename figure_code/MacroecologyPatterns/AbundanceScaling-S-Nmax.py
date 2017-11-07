from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import log10, log2, sqrt, exp, log
import scipy.optimize as opt
from math import erf, pi
import os
import sys
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def alpha2(a, N, Nmax, Nmin=1):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a))
    y += erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N
    return y

def s2(a, Nmax, Nmin=1):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10

def getNmax(N, b, slope):
    return 10 ** (b + slope*(log10(N)))

def expS(N, b, slope):
    return 10 ** (b + slope*(log10(N))) # 0.78 + 0.37*

def getS(Nrange, sb, sz, db, dz, guess, NmaxRange = [], predictNmax=True):

    Dlist = []
    Slist_ln = []
    Slist_SvN = []
    Nlist = []

    for i in range(1000):
        N = float(np.random.uniform(Nrange)[1])
        Nlist.append(N)

        Nmax = 0
        if predictNmax == True: Nmax = getNmax(N, db, dz)
        else: Nmax = np.random.uniform(NmaxRange)[1]

        Dlist.append(Nmax)
        Nmin = 1
        a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]

        S2 = s2(a, Nmax, 1)
        Slist_ln.append(S2)

        S = expS(N, sb, sz)
        Slist_SvN.append(S)

    return [log10(Slist_ln), log10(Slist_SvN), log10(Dlist), log10(Nlist)]


df = pd.read_csv(mydir + '/results/simulated_data/main/SimData.csv')
df2 = pd.DataFrame({'sim' : df['sim']})

df = df[df['species.richness'] > 0]
df = df[df['total.abundance'] > 0]
df = df[df['logmod.skew'] != 0]
df = df[df['simpson.e'] != 0]

df2['N'] = np.log10(df['total.abundance'].groupby(df['sim']).mean())
df2['D'] = np.log10(df['N.max'].groupby(df['sim']).mean())
df2['S'] = np.log10(df['species.richness'].groupby(df['sim']).mean())

df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

y1 = list(df2['D'])
y2 = list(df2['S'])
x2 = list(df2['N'])

m, b, r, p, std_err = stats.linregress(x2, y1)
db = b
dz = m

m, b, r, p, std_err = stats.linregress(x2, y2)
sb = b
sz = m

## Begin making figure
fig = plt.figure()
fs = 6
CI = 0.0001
c = '0.3'

# Adding in derived/inferred points
GO = [3.6*(10**28), 10.1*(10**28)] # estimated open ocean bacteria; Whitman et al. 1998
Pm = [2.8*(10**27), 3.0*(10**27)] # estimated Prochlorococcus; Flombaum et al. 2013
Syn = [6.7*(10**26), 7.3*(10**26)] # estimated Synechococcus; Flombaum et al. 2013

Earth = [9.2*(10**29), 31.7*(10**29)] # estimated bacteria on Earth; Kallmeyer et al. 2012
SAR11 = [2.0*(10**28), 2.0*(10**28)] # estimated percent abundance of SAR11; Morris et al. (2002)

HGx = [0.5*(10**14), 1.5*(10**14)] # estimated bacteria in Human gut; Berg (1996)
HGy = [0.05*min(HGx), 0.15*max(HGx)] # estimated most abundant bacteria in Human gut; Turnbaugh et al. (2009), & Dethlefsen et al. (2008)

COWx = [0.5*2.226*(10**15), 1.5*2.226*(10**15)] # estimated bacteria in Cow rumen; LOW:   HIGH: Whitman et al. (1998)
COWy = [0.09*min(COWx), 0.15*max(COWx)] # estimated dominance in Cow rumen; Stevenson and Weimer (2006)

Ns = []
Ss = []
DomSs = []
Ds = []

# Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
guess = 0.1019
yrange = [min(Syn), max(Pm)]
Slist_ln, Slist_SvN, Dlist, Nlist = getS(GO, sb, sz, db, dz, guess, yrange, predictNmax=False)
S_ln = np.mean(Slist_ln)
S1 = float(S_ln)
Nmax = np.mean(Dlist)
Ds.append(Nmax)
avgN = np.mean(Nlist)
Ss.append(S_ln)
GOceanS2 = float(S_ln)
N = float(avgN)
Ns.append(N)
DomSs.append(GOceanS2)

# Earth, i.e., Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
guess = 0.1060
yrange = [min(Pm), max(SAR11)]
Slist_ln, Slist_SvN, Dlist, Nlist = getS(Earth, sb, sz, db, dz, guess, yrange, predictNmax=False)
S_ln = np.mean(Slist_ln)
S1 = float(S_ln)
Nmax = np.mean(Dlist)
Ds.append(Nmax)
avgN = np.mean(Nlist)
Ss.append(S_ln)
EarthS2 = float(S_ln)
N = float(avgN)
Ns.append(N)
DomSs.append(EarthS2)

# Human Gut
guess = 0.1509
Slist_ln, Slist_SvN, Dlist, Nlist = getS(HGx, sb, sz, db, dz, guess, HGy, predictNmax=False)
S_ln = np.mean(Slist_ln)
Nmax = np.mean(Dlist)
Ds.append(Nmax)
avgN = np.mean(Nlist)
Ss.append(S_ln)
HGutS2 = float(S_ln)
N = float(avgN)
Ns.append(N)
DomSs.append(HGutS2)

# Cow Rumen
guess = 0.1
Slist_ln, Slist_SvN, Dlist, Nlist = getS(COWx, sb, sz, db, dz, guess, COWy, predictNmax=False)
S_ln = np.mean(Slist_ln)
Nmax = np.mean(Dlist)
Ds.append(Nmax)
avgN = np.mean(Nlist)
Ss.append(S_ln)
CRumenS2 = float(S_ln)
N = float(avgN)
Ns.append(N)
DomSs.append(CRumenS2)



ax = fig.add_subplot(2, 2, 1)

Micy = (np.log10(0.38)+np.arange(32)*0.93).tolist()
Micx = range(32)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = '$log$'+r'$_{10}$'+'($N_{max}$)'

d = pd.DataFrame({'N': list(x2)})
d['D'] = list(y1)
f = smf.ols('D ~ N', d).fit()

# code for prediction intervals
X0 = np.linspace(0, 5, 15)
Y0 = f.predict(exog=dict(N=X0))
X1 = np.linspace(0, 32, 100)
Y1 = f.predict(exog=dict(N=X1))
Nlist2 = X0.tolist() + x2 + X1.tolist()
Dlist2 = Y0.tolist() + y1 + Y1.tolist()

d = pd.DataFrame({'N': list(Nlist2)})
d['y'] = list(Dlist2)
f = smf.ols('y ~ N', d).fit()

st, data, ss2 = summary_table(f, alpha=CI)
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='0.7', lw=0.5, alpha=0.2)
z = np.polyfit(Nlist2, Dlist2, 1)
p = np.poly1d(z)
xp = np.linspace(0, 32, 1000)

lab = r'$S$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'

plt.xlabel(xlab, fontsize=fs+2)
plt.ylabel(ylab, fontsize=fs+2)
plt.tick_params(axis='both', labelsize=fs)
plt.xlim(0, 31)
plt.ylim(0.3, 30)

label1 = 'Abundances from published sources'
label2 = 'Simulation data: $S$ = '+str(round(10**b,2))+'$N^{'+str(round(m,2))+'}$'
label3 = 'Microbial data: $N_{max}$ = 0.38'+'$N^{0.93}$'


l1 = plt.scatter(Ns, Ds, color = '0.4', s = 10, linewidths=0.5, edgecolor='k', label=label1)
l2, = plt.plot(Micx, Micy, '--', lw=0.5, color='r', label=label3)
l3, = plt.plot(xp, p(xp), '--', lw=1, color='0.2', label=label2)

# Create a legend for the first line.
first_legend = plt.legend(handles=[l1, l2, l3], loc=2, fontsize=fs-1, frameon=False)
# Add the first legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

tdf = pd.read_csv(mydir + '/2017_09_28_1707_moretaxa.csv')
tdf = tdf[tdf['Log10Nmax'] > 0]

plot_lines = []
newL = list(tdf['Taxon'])
colors = ['purple', 'green', 'blue', 'orange', 'red']
newN = list(tdf['Log10N'])
newD = list(tdf['Log10Nmax'])

l4 = plt.scatter(newN[0], newD[0], s=10, color = colors[0], label=newL[0])
l5 = plt.scatter(newN[1], newD[1], s=30, color = colors[1], label=newL[1])
l6 = plt.scatter(newN[2], newD[2], s=10, color = colors[2], label=newL[2])
l7 = plt.scatter(newN[3], newD[3], s=10, color = colors[3], label=newL[3])
l77 = plt.scatter(newN[4], newD[4], s=10, color = colors[4], label=newL[4])

# Create another legend for the second line.
plt.legend(handles=[l4, l5, l6, l7, l77], loc=4, frameon=False, fontsize=fs-1)


ax = fig.add_subplot(2, 2, 2)

xlab = '$log$'+r'$_{10}$'+'($N$)'
ylab = '$log$'+r'$_{10}$'+'($S$)'

d = pd.DataFrame({'N': list(x2)})
d['S'] = list(y2)
f = smf.ols('S ~ N', d).fit()

# code for prediction intervals
X0 = np.linspace(0, 5, 15)
Y0 = f.predict(exog=dict(N=X0))
X1 = np.linspace(0, 32, 100)
Y1 = f.predict(exog=dict(N=X1))
Nlist2 = X0.tolist() + x2 + X1.tolist()
Slist2 = Y0.tolist() + y2 + Y1.tolist()

d = pd.DataFrame({'N': list(Nlist2)})
d['y'] = list(Slist2)
f = smf.ols('y ~ N', d).fit()

st, data, ss2 = summary_table(f, alpha=CI)
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='0.7', lw=0.5, alpha=0.2)
z = np.polyfit(Nlist2, Slist2, 1)
p = np.poly1d(z)
xp = np.linspace(0, 32, 1000)

lab = r'$S$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'

plt.xlabel(xlab, fontsize=fs+2)
plt.ylabel(ylab, fontsize=fs+2)
plt.tick_params(axis='both', labelsize=fs)
plt.xlim(0, 31)
plt.ylim(0.3, 14)

# Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
ax.text(28, 9.2, 'Ocean', fontsize=fs-2, color = 'k')

# Earth, i.e., Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
ax.text(28, EarthS2*1.05, 'Earth', fontsize=fs-2, color = 'k')

# Human Gut
ax.text(8.4, HGutS2, 'Human Gut', fontsize=fs-2, color = 'k')

# Cow Rumen
ax.text(9.3, CRumenS2*1.05, 'Cow Rumen', fontsize=fs-2, color = 'k')


label1 = 'Lognormal predictions of microbial $S$'
label2 = 'Microbial data: $S$ = 7.6'+'$N^{0.35}$'
label3 = 'Macrobial data: $S$ = 1.8'+'$N^{0.24}$'
label4 = 'Simulation data: $S$ = '+str(round(10**b,2))+'$N^{'+str(round(m,2))+'}$'

Micy = (np.log10(7.6)+np.arange(32)*0.35).tolist()
Micx = range(32)

Macy = (np.log10(1.77)+np.arange(32)*0.24).tolist()
Macx = range(32)

l8 = plt.scatter(Ns, Ss, color = '0.4', s = 10, linewidths=0.5, edgecolor='k', label=label1)
l9, = plt.plot(Micx, Micy, '--', lw=0.5, color='r', label=label2)
l10, = plt.plot(Macx, Macy, '--', lw=0.5, color='b', label=label3)
l11, = plt.plot(xp, p(xp), '--', lw=1, color='0.2', label=label4)

# Create a legend for the first line.
first_legend = plt.legend(handles=[l8, l9, l10], loc=2, fontsize=fs-1, frameon=False)
# Add the first legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

tdf = pd.read_csv(mydir + '/2017_09_28_1707_moretaxa.csv')

plot_lines = []
newL = tdf['Taxon']
colors = ['red', 'orange', 'yellow', 'green', 'blue', 'm', 'cyan', 'lime',
    'crimson', 'violet', 'Steelblue', 'brown', 'LightSeaGreen', 'Pink', 'Gold']
newN = tdf['Log10N']
newS = tdf['Log10S']

l11 = plt.scatter(newN[0], newS[0], s=20, edgecolor='0.8', linewidths=0.5, color = colors[0], label=newL[0])
l12 = plt.scatter(newN[1], newS[1], s=20, edgecolor='0.8', linewidths=0.5, color = colors[1], label=newL[1])
l13 = plt.scatter(newN[2], newS[2], s=20, edgecolor='0.8', linewidths=0.5, color = colors[2], label=newL[2])
l14 = plt.scatter(newN[3], newS[3], s=20, edgecolor='0.8', linewidths=0.5, color = colors[3], label=newL[3])
l15 = plt.scatter(newN[4], newS[4], s=20, edgecolor='0.8', linewidths=0.5, color = colors[4], label=newL[4])
l16 = plt.scatter(newN[5], newS[5], s=20, edgecolor='0.8', linewidths=0.5, color = colors[5], label=newL[5])
l17 = plt.scatter(newN[6], newS[6], s=20, edgecolor='0.8', linewidths=0.5, color = colors[6], label=newL[6])
l18 = plt.scatter(newN[7], newS[7], s=20, edgecolor='0.8', linewidths=0.5, color = colors[7], label=newL[7])
l19 = plt.scatter(newN[8], newS[8], s=20, edgecolor='0.8', linewidths=0.5, color = colors[8], label=newL[8])
l20 = plt.scatter(newN[9], newS[9], s=20, edgecolor='0.8', linewidths=0.5, color = colors[9], label=newL[9])
l21 = plt.scatter(newN[10], newS[10], s=20, edgecolor='0.8', linewidths=0.5, color = colors[10], label=newL[10])
l22 = plt.scatter(newN[11], newS[11], s=20, edgecolor='0.8', linewidths=0.5, color = colors[11], label=newL[11])

# Create another legend for the second line.
plt.legend(handles=[l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22],
    loc=2, frameon=False, fontsize=fs, bbox_to_anchor=(1.05, 1),
    borderaxespad=0.)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/S-AbundanceScaling3.png', dpi=600, bbox_inches = "tight")
plt.close()
