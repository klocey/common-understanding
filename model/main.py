from __future__ import division
from random import shuffle, choice, randint, sample
from os.path import expanduser
from numpy import log10
from scipy import stats
import numpy as np
import time
import copy
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/common-understanding/model")
GenPath = mydir + "GitHub/common-understanding/results/simulated_data/"

col_headers = 'sim,ct,num-firms,info-lag,commons-growth,commons-size,firm-total-size,firm-avg-size'
OUT = open(GenPath + 'SimData-main.csv', 'w+')
print>>OUT, col_headers
OUT.close()


def output(iD, sD, rL, sim, ct, ps):
    r, u, gr, mt, q = ps
    IndIDs, SpIDs = [], []
    for k, v in iD.items():
            IndIDs.append(k)
            SpIDs.append(v['sp'])

    N = len(IndIDs)
    R = sum(rL)
    S = len(list(set(SpIDs)))

    if N > 0:
        OUT = open(GenPath + 'SimData-MainR.csv', 'a')
        outlist = [sim, u, r, gr, mt, q, ct, N, S, ES, Nm, lms]
        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

    print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  N, '  S:', '%4s' %  S,  '  R:', '%4s' %  R
    return


def firms(fD, ps, num_firms, num_res):
    take_max, give_max, res_in_max = ps

    p = time.time()
    fD[p] = {'take' : np.random.uniform(take_max, 0)}
    fD[p]['give'] = np.random.uniform(give_max, 0)
    fD[p]['q'] = 0
    fD[p]['res'] = choice(range(0, num_res))
    return fD


def commons(rD, ps, num_res):
    take_max, give_max, res_in_max = ps
    for i in range(num_res):
        rD[i] = {'res' : i}
        rD[i]['val'] = np.random.uniform(10, 1000)
    return rD


def commons_growth(rD, ps):
    take_max, give_max, res_in_max = ps
    for k, v in rD.items():
        rD[k] += np.random.uniform(0, res_in_max)
    return rD


def take(fD, rD, ps):
    r, u, gr, mt, q = ps
    keys = list(fD)
    shuffle(keys)

    for k in keys:
        r = fD[k]['res']
        if rD[r]['val'] > 0:
            fD[k]['q'] += min([rD[r]['val'], 1])
            rD[r] -= min([rD[r]['val'], 1])
    return fD, rD


def give(fD, rD, ps):
    r, u, gr, mt, q = ps
    keys = list(rD)
    shuffle(keys)

    for k in keys:
        r = rD[k]['res']
        if fD[r]['val'] > 0:
            rD[k]['q'] += min([fD[r]['val'], 1])
            fD[r] -= min([fD[r]['val'], 1])
    return fD, rD


def iter_procs(fD, rD, ps, ct):
    procs = range(3)
    shuffle(procs)
    for p in procs:
        if p == 0: rD = commons_growth(rD, ps)
        elif p == 1: fD, rD = take(fD, rD, ps)
        elif p == 2: fD, rD = give(fD, rD, ps)
    return [fD, rD, ct+1]


for sim in range(10**5):
    print '\n'
    take_max = 10**np.random.uniform(-2, -1)
    give_max = 10**np.random.uniform(-2, -1)
    res_in_max = 10**np.random.uniform(-1, 1)
    num_res = choice(range(1,10))
    num_firms = choice(range(1,10))

    fD, rD, ct = {}, {}, 0
    ps = take_max, give_max, res_in_max
    fD = firms(fD, rD, ps, num_firms, num_res)
    rD = commons(rD, ps, num_res)

    while ct < 300:
        fD, rD, ct = iter_procs(fD, rD, ps, ct)
        if ct > 200 and ct%10 == 0: output(fD, rD, sim, ct, ps)
