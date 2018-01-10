#! /usr/bin/env python

import numpy as np
import sys
import fitters
import gvar as gv



def bs_corr_list(corr,Nbs,Mbs = None,seed=None,return_b0 = True):
    #feed in correlator data in format (N_cfgs, ...)
    #Andre wrote this
    #this returns an array of Nbs lists of Mbs random indices in the correlator data
    #you can use this to generate a new bootstrapped data set from the original data
    if Mbs is None:
        Mbs = corr.shape[0]
    else:
        print "Currently only handles Mbs = N_cfgs"
        sys.exit(-1)
    np.random.seed(seed) # if None - it does not seed - I checked 14 May 2013
    # make bs_lst of shape (Nbs,Mbs)
    bs_list = np.random.randint(0,corr.shape[0],(Nbs,Mbs))
    original = np.arange(corr.shape[0])
    #print original
    if return_b0:
        bs_list = np.insert(bs_list,0,original,axis=0)
        #print "return b0"
    return bs_list
import numpy as np
import sys

def bs_corr_list(corr,Nbs,Mbs = None,seed=None,return_b0 = True):
    #feed in correlator data in format (N_cfgs, ...)
    #this returns an array of Nbs lists of Mbs random indices in the correlator data
    #you can use this to generate a new bootstrapped data set from the original data
    if Mbs is None:
        Mbs = corr.shape[0]
    else:
        print "Currently only handles Mbs = N_cfgs"
        sys.exit(-1)
    np.random.seed(seed) # if None - it does not seed - I checked 14 May 2013
    # make bs_lst of shape (Nbs,Mbs)
    bs_list = np.random.randint(0,corr.shape[0],(Nbs,Mbs))
    original = np.arange(corr.shape[0])
    #print original
    if return_b0:
        bs_list = np.insert(bs_list,0,original,axis=0)
        #print "return b0"
    return bs_list

def s_func(psq,L,Nlam=19):
    '''
    The S-function is equivalent to the Riemann-Zeta function in the A1+ irrep with P_com=0
    I assume that psq is a number or a 1D np array or list
    '''
    sqrsTbl = [6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24]
    x = psq * L**2 / 4 / np.pi**2
    try:
        Nbs = len(psq)
    except:
        Nbs = 1
    if Nbs == 1:
        S = sum([x**7 * sqrsTbl[i]/((i+1)**7 * (i+1-x)) for i in range(int(Nlam)+1)])
    else:
        S = np.zeros([Nbs])
        for bs in range(Nbs):
            S[bs] = sum([x[bs]**7 * sqrsTbl[i]/((i+1)**7 * (i+1-x[bs])) for i in range(Nlam)])
    S -= 1/x
    S -= 2.8373*np.pi
    S += 16.5323*x
    S += 8.40192 * x**2
    S += 6.94581 * x**3
    S += 6.4261191 * x**4
    S += 6.20215 * x**5
    S += 6.09818 * x**6
    S = S / np.pi / L
    return(S)

def psq_int(E_2pi,m):
    '''
    compute psq from 2 particle energer and mass
    E_2pi = 2 * sqrt( m**2 + p**2)
    psq = E_2pi**2 / 4 - m**2
    '''
    return 0.25 * E_2pi**2 - m**2

def pcotd(E_2pi,m,L):
    '''
    compute p cot(delta) in lattice units given the 2-particle energy and mass
    '''
    psq   = psq_int(E_2pi,m)
    pcotd = s_func(psq,L)
    return pcotd
 