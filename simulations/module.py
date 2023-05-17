#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 12:52:01 2021

@author: msanter
"""

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.stats import hypergeom,betabinom
import pickle
import collections.abc

# initial population with a fraction of f cell with one mutant copy
def Ninitfun(n, f, Npop):
    Ninit = np.zeros(n + 1)
    Ninit[0] = 1 - f
    Ninit[1] += f
    Ninit *= Npop
    return Ninit

def betabinom_pdf_special(k,n,a,b):
    if a==0 and b==0:
        raise Exception("a=b=0")
    elif a==0:
        return float(k==0)
    elif b==0:
        return float(k==n)
    else:
        return betabinom.pmf(k,n,a,b)

def stochbottleSim(
    n,  # replicon copy number
    rmt=1.,  # mutant growth rate
    rwt=1. - 0.1,  # wild-type growth rate
    r_array=None,
    Td=0,  # time per day, if Td=0 the equilibrium Td->inf is calculated
    Ninit=None, # initial composition of cells, 
                # if Ninit==None : initial population 
                # has frequency f of mutant cells with one mt. copy (rest wt.)
    Npop=1e7, # (if Ninit==None) initial population size
    f=1e-4, # (if Ninit==None) mutant cells with 1 mutant copy
    rep="reg",  # replication/segregation mode
    D=24,  # days in the experiment
    Nc=1e9,  # carrying capacity (scalar or array possible)
    b=0.01,  # bottleneck factor
    file='', # file for saving output
    eqthr=.01, # thresfold for equilibrium (is set below 1 to get good equilibrium)
    seed=0, # seed value for stochastic simulation
    d0=1, # starting day of simulation
    swb=False, # True: starts with bottleneck (if True Nc[0]*b is used as bottleneck
                    # and Nc[1:] used for carrying capacities )
    p=None, # plasmid segregation matrix
    pget=False, # return segregation matrix only
    verbose=False, # print population states during simulation
    establishment=False, # if True, stop simulation if establishment is reached 
):

    np.random.seed(seed)

    # convert Nc to an array of constant carrying capacity for all days
    
    if not (isinstance(Nc, collections.abc.Sequence) or isinstance(Nc, np.ndarray)):
        Nc=np.ones(D+1)*Nc

    # initialize population array for the first day (start of day)
    N = Ninitfun(n,f,Npop) if Ninit==None else Ninit

    d=d0 # day

    # timeseries of population compositions at the start of days after bottleneck
    ts_sod = []  
    # timeseries of population compositions at the end of the days after growth
    ts_eod = []
    # list of intra-day population dynamics
    ts_intraday = [] 

    # initialize stochastic matrix p of replicon inheritance
    # (expected number of j-type offspring at i-type cell division)
    if not(p is None) and pget==False:
        pass
    elif rep=='reg':
    # implementation of regular replication and random segregation
        p = np.array(
            [
                [
                    2 * hypergeom(2 * n, 2 * i, n).pmf(j)
                    for j in range(0, n + 1)
                ]
                for i in range(0, n + 1)
            ]
        )
    elif rep=='ran':
    # implementation of random replication and random segregation
        def pranrep(i,j):
            tmp=0
            for k in range(i,i+n+1):
                tmp+= betabinom_pdf_special(k - i, n, i, n - i) \
                    * 2 * hypergeom(2 * n, k, n).pmf(j)
            return tmp
        p = np.array(
            [
                [
                    pranrep(i,j)
                    for j in range(0, n + 1)
                ]
                for i in range(0, n + 1)
            ]
        )

    # cell-type growth function r(i)
    # implementation of a dominant mutation
    def r(i):
        if r_array==None:
            if i > 0:
                return rmt
            else:
                return rwt
        else:
            return r_array[i]

    # define function to evaluate if the equilibrium is reached

    # inter-day population dynamics (stochastic)
    while True:
        if d>D:
            break
        if d!=d0 or swb:
            # apply bottleneck 
            pvals = N / sum(N)
            # Avoid pvals to get negative by division of 0 by large number
            pvals[np.argwhere(pvals<0)]=0 
            N = np.random.multinomial(round(Nc[d-d0] * b), pvals)
        ts_sod.append(N.copy()) # store pop. state
        # define intra-day population dynamics 

        def f(t, N):
            dNdt = np.zeros(n + 1)
            for i in range(n + 1):
                for j in range(n + 1):
                    a=r(j) * (1 - np.sum(N) / Nc[d-d0+swb]) * (p[j, i] - (i == j)) * N[j]
                    dNdt[i] += a
            return dNdt
        # define function for determining the equilibrium
        def hit_eq(t, X):
            # ignore division by 0 in some entries
            # maximum of relative derivatives
            delta=np.max(np.abs(f(t,X)))
            return delta-eqthr
        hit_eq.terminal = True
        # compute the population growth dynamics for the current day

        sol=None
        if Td!=0:
            sol = solve_ivp(f, [0, Td], N, method='BDF', dense_output=True)
            N = sol.y.T[-1]  
        # shortcut for homozygous populations
        elif (N[0]==np.sum(N[:])):
            N[0]=Nc[d-d0+swb]
            N[1:]=0
        elif (N[-1]==np.sum(N[:])):
            N[-1]=Nc[d-d0+swb]
            N[:-1]=0
        else: # integration
            sol = solve_ivp(f, [0, 1e10], N, method='Radau', dense_output=True,
                events=[hit_eq] )
            if len(sol.t_events[0])==0:
                raise ValueError('Warning: Equilibrium not reached')
            # population at the end of the day
            N = sol.y.T[-1]  
        # If heterozygotes frequency was zero at the start of the day then set
        # abundance of all heterozygous types at the end of the day as 
        # integration results in negative cell-type abundances in rare cases.
        if np.sum(N[1:n])==0:
            for i in range(1,n):
                N[i]=0
        # copy population state to list of Time Series of intra-day dynamics
        if sol!=None:
            ts_intraday.append([sol.t,sol.y.T])
        else:
            ts_intraday.append(sol)
        # copy population state to list of Time Series at the End Of the Day
        ts_eod.append(N)

        if verbose:
            print(d,ts_sod[-1],np.sum(ts_sod[-1]),ts_eod[-1],np.sum(ts_eod[-1]))

        if establishment:
            if (1-N[-1]/Nc[d-d0])**(b*Nc[d-d0])<=.0001: # probability of stochastic loss 
                                            # at bottleneck is tiny ->establishment
                break
            if N[0]==sum(N): # extinction of the mutant
                break

        d+=1

    if file!='':
            pickle.dump((np.array(ts_sod), np.array(ts_eod), ts_intraday),
            open(file,'wb+'))
    if pget==True:
        return p
    else:
        return np.array(ts_sod), np.array(ts_eod), ts_intraday