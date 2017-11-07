from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
import scipy as sc
from scipy import stats

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
GenPath = mydir + "GitHub/ScaleEmerge/results/simulated_data/"

col_headers = 'frac', 'choose', 'Eslope', 'Eint', 'Dslope', 'Dint', 'Rslope', 'Rint', 'Sslope', 'Sint'
OUT = open(mydir + 'SimData-Random-ModelResults.csv', 'w+')
print>>OUT, col_headers
OUT.close()

df = pd.read_csv(mydir + '/results/simulated_data/SimData-Random.csv')
df = df[df['species.richness'] > 0]
df = df[df['total.abundance'] > 0]
df = df[df['logmod.skew'] != 0]
df = df[df['simpson.e'] != 0]
df = df[df['sim'] < 10]

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
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()


dfo = pd.DataFrame(columns=())

clrs = list(df2['color'])
unique_colors = list(set(clrs))

for i, clr in enumerate(unique_colors):
    df4 = pd.DataFrame({'clr' : df2['color']})
    df4['N'] = df2['N']
    df4['S'] = df2['S']
    df4['E'] = df2['E']
    df4['R'] = df2['R']
    df4['D'] = df2['D']
    df4['frac'] = df2['frac']
    df4['choose'] = df2['choose']
    df4 = df4[df4['clr'] == clr]

    outlist = [df4['frac'][0], df4['choose'][0]]
    Em, Eb, r, p, std_err = stats.linregress(df4['N'], df4['E'])
    Dm, Db, r, p, std_err = stats.linregress(df4['N'], df4['D'])
    Rm, Rb, r, p, std_err = stats.linregress(df4['N'], df4['R'])
    Sm, Sb, r, p, std_err = stats.linregress(df4['N'], df4['S'])
    outlist.extend([Em, Eb, Dm, Db, Rm, Rb, Sm, Sb])

    OUT = open(GenPath + 'SimData-Random-ModelResults.csv', 'a')
    outlist = str(outlist).strip('[]')
    outlist = outlist.replace(" ", "")
    print>>OUT, outlist
    OUT.close()
