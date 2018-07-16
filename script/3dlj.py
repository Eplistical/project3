#!/usr/bin/env python3

# repeat Sec. II.A results in van Swol's paper 10.1063/1.464945
# 3d LJ bulk fluid w/ GCMC & GCMD methods

import numpy as np
from numpy import linalg as la
from numba import jit
import argparse
from matplotlib import pyplot as plt


# config
kT = 1.0
beta = 1.0 / kT
mu = -3.2
V = 512
L = V**(1.0 / 3.0)
Vc = V
Lc = Vc**(1.0 / 3.0)
dxmax = 1.0

Nrun = 1
Nstep_eql = int(1e5)
Nstep_run = int(3e5)


class LJ_with_cutoff(object):
    """pair-wise potential"""
    rc = 2.5
    U_rc = None

    def __init__(self):
        self.U_rc = self.LJ_raw(self.rc)

    def __call__(self, r):
        if r <= self.rc:
            return self.LJ_raw(r) - self.U_rc
        else:
            return 0.0

    #@jit
    def LJ_raw(self, r):
        assert r > 0.0
        tmp = 1.0 / r**6
        return 4 * tmp * (tmp - 1)

LJ = LJ_with_cutoff()


def parse_args():
    global Nstep_eql, Nstep_run, Nrun
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nrun', type=int, help="run times")
    parser.add_argument('--Nstep_eql', type=int, help="Nstep for equilibrium")
    parser.add_argument('--Nstep_run', type=int, help="Nstep for running")
    args = parser.parse_args()
    if args.Nstep_eql is not None:
        Nstep_eql = args.Nstep_eql
    if args.Nstep_run is not None:
        Nstep_run = args.Nstep_run
    if args.Nrun is not None:
        Nrun = args.Nrun


parse_args()
print('# V = %.2f, Vc = %.2f' %(V, Vc))
print('# Nstep_eql = %d, Nstep_run = %d, Nrun = %d' % (Nstep_eql, Nstep_run, Nrun))


# functions
#@jit
def cal_U(swarm):
    N = cal_N(swarm)
    rst = 0.0
    for i in range(N - 1):
        for j in range(i + 1, N):
            rst += LJ(la.norm(swarm[i] - swarm[j]))
    return rst


#@jit
def cal_N(swarm):
    return len(swarm)


#@jit
def is_in_Vc(ptcl):
    return np.all(np.logical_and(ptcl < Lc, ptcl > 0.0))


#@jit
def cal_Nc(swarm):
    rst = 0
    for ptcl in swarm:
        if is_in_Vc(ptcl):
            rst += 1
    return rst


#@jit
def rand_in_Vc():
    """generate a rancom position in Lc"""
    return np.random.rand(3) * Lc


#@jit
def MC_step(swarm):
    swarm = shuffle(swarm)
    if np.random.rand() < 0.5:
        swarm = create(swarm)
    else:
        swarm = destruct(swarm)
    return swarm


#@jit
def shuffle(swarm):
    new_swarm = list()
    for k, ptcl in enumerate(swarm):
        dx = (np.random.rand(3) * 2.0 - 1.0) * dxmax
        new_ptcl = np.copy(ptcl) + dx
        # periodic condition
        for i in range(3):
            if new_ptcl[i] < 0.0:
                new_ptcl[i] += L
            elif new_ptcl[i] > L:
                new_ptcl[i] -= L
        assert np.all(np.logical_and(new_ptcl >= 0.0, new_ptcl <= L))
        # calc dU
        dU = 0.0
        for j, ptcl_j in enumerate(swarm):
            if j != k:
                dU += LJ(la.norm(new_ptcl - ptcl_j)) - LJ(la.norm(ptcl - ptcl_j))
        # make decision
        prob = np.exp(-beta * dU)
        if np.random.rand() < prob:
            new_swarm.append(new_ptcl)
        else:
            new_swarm.append(ptcl)
    return new_swarm


#@jit
def create(swarm):
    new_ptcl = rand_in_Vc()
    # cal dU
    dU = 0.0
    for i, ptcl in enumerate(swarm):
        dU += LJ(la.norm(new_ptcl - ptcl))
    # make decision
    Nc = cal_Nc(swarm)
    prob = np.exp(-beta * dU + beta * mu) * Vc / (Nc + 1)
    if np.random.rand() < prob:
        swarm.append(new_ptcl)
    return swarm


#@jit
def destruct(swarm):
    Nc = cal_Nc(swarm)
    if Nc > 0:
        # get idx of the ptcl to destruct
        tmp = np.random.choice(Nc)
        idx = 0
        while tmp > 0:
            if is_in_Vc(swarm[idx]):
                tmp -= 1
            idx += 1
        # calc dU
        dU = 0.0
        for i, ptcl in enumerate(swarm):
            if i != idx:
                dU -= LJ(la.norm(ptcl - swarm[idx]))
        # make decision
        prob = np.exp(-beta * dU - beta * mu) / Nc * Vc
        if np.random.rand() < prob:
            del swarm[idx]
    return swarm


#@jit
def main():
    for irun in range(Nrun):
        swarm = list()
        swarm.append(np.array([L/2, L/2, L/2]))

        # equilibrate
        for istep in range(Nstep_eql):
            print(istep)
            swarm = MC_step(swarm)

        # collect data
        for istep in range(Nstep_run):
            print(istep)
            swarm = MC_step(swarm)

    # output



if __name__ == '__main__':
    main()


# END
