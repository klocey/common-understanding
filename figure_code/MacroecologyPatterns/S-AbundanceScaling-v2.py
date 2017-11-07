from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import log10, log2, sqrt, exp, log
import scipy.optimize as opt
from math import erf, pi
import os
import sys
import linecache
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table


def alpha2(a, N, Nmax, Nmin=1):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a))
    y += erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N
    return y

def s2(a, Nmax, Nmin=1):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10 of Curtis and Sloan (2002)

def getNmax(N, b, slope):
    return 10 ** (b + slope*(log10(N)))

def expS(N, b, slope):
    return 10 ** (b + slope*(log10(N)))

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


c = '0.3'
datasets = []
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
tail = '-SADMetricData.txt'
GoodNames = ['EMPclosed', 'HMP', 'BIGN', 'TARA', 'BOVINE', 'HUMAN', 'LAUB', 'SED', 'CHU',
    'CHINA', 'CATLIN', 'FUNGI', 'HYDRO', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA']

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
for name in os.listdir(mydir +'data/micro'):
    if name in GoodNames: pass
    else: continue
    path = mydir+'data/micro/'+name+'/'+name+tail
    num_lines = sum(1 for line in open(path))
    datasets.append([name, 'micro', num_lines])
    print name, num_lines

for name in os.listdir(mydir +'data/macro'):
    if name in GoodNames: pass
    else: continue
    path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
    num_lines = sum(1 for line in open(path))
    datasets.append([name, 'macro', num_lines])
    print name, num_lines

fs = 12 # font size used across figures
MicIntList, MicCoefList, MacIntList, MacCoefList, R2List, metlist = [[], [], [], [], [], []]
Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

for dataset in datasets:
    lines = []
    name, kind, numlines = dataset
    small = ['BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED']
    big = ['HUMAN', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO']

    if kind == 'macro':
        lines = np.random.choice(range(1, numlines+1), 1000, replace=True)
    elif name in small:
        lines = np.random.choice(range(1, numlines+1), 200, replace=True)
    elif name in big:
        lines = np.random.choice(range(1, numlines+1), 500, replace=True)
    elif name == 'TARA':
        lines = np.random.choice(range(1, numlines+1), 500, replace=True)
    else:
        lines = np.random.choice(range(1, numlines+1), 500, replace=True)
    if kind == 'micro': path = mydir+'data/'+kind+'/'+name+'/'+name+tail
    else: path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

    for line in lines:
        data = linecache.getline(path, line)
        radDATA.append(data)


for data in radDATA:
    data = data.split()
    name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data
    Nlist.append(np.log10(float(N)))
    Slist.append(np.log10(float(S)))
    KindList.append(kind)
    if kind == 'micro': klist.append('b')
    if kind == 'macro': klist.append('r')


# Multiple regression
d = pd.DataFrame({'N': list(Nlist)})
d['S'] = list(Slist)
d['Kind'] = list(KindList)
f = smf.ols('S ~ N * Kind', d).fit()

# Microbes
X1 = np.linspace(0, 32, 1000)
mKind = ['micro']*1000
Y1 = f.predict(exog=dict(N=X1, Kind=mKind))
Mic_Nlist2 = X1.tolist()
Mic_Slist2 = Y1.tolist()

d = pd.DataFrame({'N': list(Mic_Nlist2)})
d['S'] = list(Mic_Slist2)
mf = smf.ols('S ~ N', d).fit()

st, data, ss2 = summary_table(mf, alpha=0.05)
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T
#plt.fill_between(Mic_Nlist2, pred_ci_low, pred_ci_upp, color='Skyblue', lw=0.0, alpha=0.6)
plt.plot(Mic_Nlist2, Mic_Slist2,  color='c', ls='--', lw=0.5, alpha=1)

# Macrobes
X1 = np.linspace(0, 32, 1000)
mKind = ['macro']*1000
Y1 = f.predict(exog=dict(N=X1, Kind=mKind))
Mac_Nlist2 = X1.tolist()
Mac_Slist2 = Y1.tolist()

d = pd.DataFrame({'N': list(Mac_Nlist2)})
d['S'] = list(Mac_Slist2)
mf = smf.ols('S ~ N', d).fit()

st, data, ss2 = summary_table(mf, alpha=0.05)
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T
#plt.fill_between(Mac_Nlist2, pred_ci_low, pred_ci_upp, color='LightCoral', lw=0.0, alpha=0.6)
plt.plot(Mac_Nlist2, Mac_Slist2,  color='m', ls=':', lw=0.5, alpha=1)


MacIntList.append(f.params[0])
MacCoefList.append(f.params[2])
if f.pvalues[1] < 0.05: MicIntList.append(f.params[1] + f.params[0])
else: MicIntList.append(f.params[0])

if f.pvalues[3] < 0.05: MicCoefList.append(f.params[3] + f.params[2])
else: MicCoefList.append(f.params[2])

MicInt = round(np.mean(MicIntList), 2)
MicCoef = round(np.mean(MicCoefList), 2)
MacInt = round(np.mean(MacIntList), 2)
MacCoef = round(np.mean(MacCoefList), 2)



mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')

df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
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

dz, db, r, p, std_err = stats.linregress(x2, y1)
sz, sb, r, p, std_err = stats.linregress(x2, y2)

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

st, data, ss2 = summary_table(f, alpha=0.05)
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='0.7', lw=0.5, alpha=0.2)
z = np.polyfit(Nlist2, Slist2, 1)
p = np.poly1d(z)
xp = np.linspace(0, 32, 1000)

plt.xlabel(xlab, fontsize=fs+2)
plt.ylabel(ylab, fontsize=fs+2)
plt.tick_params(axis='both', labelsize=fs)
plt.xlim(0, 31)
plt.ylim(0.3, 14)


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

# Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
guess = 0.1019
yrange = [min(Syn), max(Pm)]
Slist_ln, Slist_SvN, Dlist, Nlist = getS(GO, sb, sz, db, dz, guess, yrange, predictNmax=False)
S_ln = np.mean(Slist_ln)
S1 = float(S_ln)
Nmax = np.mean(Dlist)
avgN = np.mean(Nlist)
Ss.append(S_ln)
S2 = float(S_ln)
N = float(avgN)
ax.text(28, 9.2, 'Ocean', fontsize=fs-4, color = 'k')
Ns.append(N)
DomSs.append(S2)

# Earth, i.e., Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
guess = 0.1060
yrange = [min(Pm), max(SAR11)]
Slist_ln, Slist_SvN, Dlist, Nlist = getS(Earth, sb, sz, db, dz, guess, yrange, predictNmax=False)
S_ln = np.mean(Slist_ln)
S1 = float(S_ln)
Nmax = np.mean(Dlist)
avgN = np.mean(Nlist)
Ss.append(S_ln)
S2 = float(S_ln)
N = float(avgN)
ax.text(28, S2*1.05, 'Earth', fontsize=fs-4, color = 'k')
Ns.append(N)
DomSs.append(S2)

# Human Gut
guess = 0.1509
Slist_ln, Slist_SvN, Dlist, Nlist = getS(HGx, sb, sz, db, dz, guess, HGy, predictNmax=False)
S_ln = np.mean(Slist_ln)
Nmax = np.mean(Dlist)
avgN = np.mean(Nlist)
Ss.append(S_ln)
S2 = float(S_ln)
N = float(avgN)
ax.text(7, S2, 'Human Gut', fontsize=fs-4, color = 'k')
Ns.append(N)
DomSs.append(S2)

# Cow Rumen
guess = 0.1
Slist_ln, Slist_SvN, Dlist, Nlist = getS(COWx, sb, sz, db, dz, guess, COWy, predictNmax=False)
S_ln = np.mean(Slist_ln)
Nmax = np.mean(Dlist)
avgN = np.mean(Nlist)
Ss.append(S_ln)
S2 = float(S_ln)
N = float(avgN)
ax.text(9, S2*1.05, 'Cow Rumen', fontsize=fs-4, color = 'k')
Ns.append(N)
DomSs.append(S2)

label1 = 'Lognormal predictions of microbial $S$'
label2 = 'Simulation data: $S$ = '+str(round(10**sb,2))+'$N^{'+str(round(sz,2))+'}$'
label3 = 'Microbial data: $S$ = '+str(round(10**MicInt,2))+'$N^{'+str(round(MicCoef,2))+'}$'
label4 = 'Macrobial data: $S$ = '+str(round(10**MacInt,2))+'$N^{'+str(round(MacCoef,2))+'}$'

Micy = (MicInt+np.arange(32)*MicCoef).tolist()
Micx = range(32)
Macy = (MacInt+np.arange(32)*MacCoef).tolist()
Macx = range(32)

l1 = plt.scatter(Ns, Ss, color = '0.4', s = 50, linewidths=2, edgecolor='k', label=label1)
l2, = plt.plot(xp, p(xp), '--', lw=1, color='0.3', label=label2)
l3, = plt.plot(Micx, Micy, '--', lw=0.5, color='b', label=label3)
l4, = plt.plot(Macx, Macy, '--', lw=0.5, color='r', label=label4)
#plt.hexbin(x2, y2, mincnt=1, gridsize = 20, bins='log', cmap=plt.cm.jet)

# Create a legend for the first line.
first_legend = plt.legend(handles=[l1, l2, l3, l4], loc=2, frameon=False)
# Add the first legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

tdf = pd.read_csv(mydir + '/2017_09_28_1707_moretaxa.csv')

plot_lines = []
newL = tdf['Taxon']
colors = ['red', 'orange', 'yellow', 'green', 'blue', 'm', 'cyan', 'lime','crimson','violet','Steelblue']
newN = tdf['Log10N']
newS = tdf['Log10S']

l4 = plt.scatter(newN[0], newS[0], s=50, color = colors[0], label=newL[0])
l5 = plt.scatter(newN[1], newS[1], s=50, color = colors[1], label=newL[1])
l6 = plt.scatter(newN[2], newS[2], s=50, color = colors[2], label=newL[2])
l7 = plt.scatter(newN[3], newS[3], s=50, color = colors[3], label=newL[3])
l8 = plt.scatter(newN[4], newS[4], s=50, color = colors[4], label=newL[4])
l9 = plt.scatter(newN[5], newS[5], s=50, color = colors[5], label=newL[5])
l10 = plt.scatter(newN[6], newS[6], s=50, color = colors[6], label=newL[6])
l11 = plt.scatter(newN[7], newS[7], s=50, color = colors[7], label=newL[7])
l12 = plt.scatter(newN[8], newS[8], s=50, color = colors[8], label=newL[8])
l13 = plt.scatter(newN[9], newS[9], s=50, color = colors[9], label=newL[9])
l14 = plt.scatter(newN[10], newS[10], s=50, color = colors[10], label=newL[10])

# Create another legend for the second line.
plt.legend(handles=[l4, l5, l6, l7, l8, l9, l10, l11,l12,l13,l14], loc=4, frameon=False, fontsize=fs-4)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/S-AbundanceScaling.png', dpi=600, bbox_inches = "tight")
plt.close()
