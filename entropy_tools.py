"""
Functions for doing mutual information analysis on expression and some predictor
"""

# all the useful packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os
import cPickle as pickle
from sys import exit
import corner
from scipy.optimize import minimize
import emcee

import seaborn as sns

rc = {'lines.linewidth': 10,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)

def Shannon(prob):
    """
    Shannon entropy associated with probability prob
    """
    if prob==1:
        return 0
    if prob==0:
        return 0
    return -(prob * np.log2(prob) + (1-prob) * np.log2(1-prob))

def p_k (params, E_k):
    """
    Our probability model
    """
    a,b,c = params
    ps =  1.0 / (c + np.exp(b*E_k + a))
    return np.array([min((1.0, p)) for p in ps])

def log_likelihood (params, Es, Fs):
    """
    The log likelihood.
    """
    a,b,c, = params
    pks = p_k (params, Es)
    pks = np.array([p if F==1 else 1-p for p,F in zip(pks, Fs)])
    return np.sum(np.log(pks))

def log_posterior (params, Es, Fs):
    """
    Use the log-likelihood, constant priors on a and c, and a
    Jeffreys prior on b to calculate the log_posterior
    """
    a,b,c = params
    if c < 0 or c > 1:
        return -np.inf
    if b < 0 or b > 5:
        return -np.inf
    if a <= -5 or a > 5:
        return -np.inf
    post = log_likelihood (params, Es, Fs) - np.log(b)
    return post

def negative_lot_post (params, Es, Fs):
    return -log_posterior (params, Es, Fs)

if __name__ == "__main__":
    with open ('2015-10-20_all_res.pkl') as f:
        df = pickle.load(f)
    n = float(len(df.index))
    fraction_expressed = (np.sum(df['expression']) + n) / n / 2.

    # mcmc to sample posterior
    n_dim = 3        # number of parameters in the model
    n_walkers = 50   # number of MCMC walkers
    n_burn = 500     # "burn-in" period to let chains stabilize
    n_steps = 50000   # number of MCMC steps to take after burn-in


    # p0[i,j] is the starting point for walk i along variable j.
    p0 = np.empty((n_walkers, n_dim))
    p0[:,0] = np.random.uniform(0, 5, n_walkers)           # a
    p0[:,1] = np.random.exponential(0.1, n_walkers)             # b
    p0[:,2] = np.random.uniform(0, 1, n_walkers)   # c
    args = (df['M'], df['expression'])
#     PS = [log_posterior(pos,df['E'], df['expression']) for pos in p0]
#     print PS
#     exit('')
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_posterior,
                                args=args)
    print 'Burn in...'
    pos, prob, state = sampler.run_mcmc(p0, n_burn, storechain=False)
    print 'Run mcmc...'
    _ = sampler.run_mcmc(pos, n_steps)
    name = '2015-10-29_information_mcmc_M'
    with open(name+'.pkl','w') as f:
        pickle.dump (sampler, f)


    fig0, ax = plt.subplots(n_dim, 1, sharex=True)
    for i in range(n_dim):
        ax[i].plot(sampler.chain[0,:,i], 'k-', lw=0.2)
        ax[i].plot([0, n_steps-1],
                   [sampler.chain[0,:,i].mean(), sampler.chain[0,:,i].mean()], 'r-')

    ax[1].set_xlabel('sample number')
    ax[0].set_ylabel('$a$')
    ax[1].set_ylabel('$b$')
    ax[2].set_ylabel('$c$')

    fig = corner.corner(sampler.flatchain,
                        labels=[r'$a$', r'$b$', r'$c$'], bins=500,plot_contours=True)
    fig.savefig(name+'_corner.pdf')
    fig0.savefig(name+'_chains.pdf')
    plt.show()


