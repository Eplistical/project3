#!/usr/bin/env python3

# repeat Sec. II.A results in van Swol's paper 10.1063/1.464945
# 1d hard rod w/ GCMC method

import numpy as np
from numba import jit
import argparse


# config
Lx = 2.0
A = 0.0
B = 2.0
beta = 1.0
mu = 1.0
dxmax = 0.5

Nstep_run = int(10e6)
Nstep_eql = int(5e5)
Nrun = 1


def parse_args():
    global A, B, Nstep_eql, Nstep_run, Nrun
    parser = argparse.ArgumentParser()
    parser.add_argument('--A', type=float, help="control volume begin")
    parser.add_argument('--B', type=float, help="control volume end")
    parser.add_argument('--Nrun', type=int, help="run times")
    parser.add_argument('--Nstep_eql', type=int, help="Nstep for equilibrium")
    parser.add_argument('--Nstep_run', type=int, help="Nstep for running")
    args = parser.parse_args()
    if args.A is not None:
        A = args.A
    if args.B is not None:
        B = args.B
    assert 0 <= A and A <= B and B <= Lx
    if args.Nstep_eql is not None:
        Nstep_eql = args.Nstep_eql
    if args.Nstep_run is not None:
        Nstep_run = args.Nstep_run
    if args.Nrun is not None:
        Nrun = args.Nrun

parse_args()
print('# Lx = %.2f, A = %.2f, B = %.2f' %(Lx, A, B))
print('# Nstep_eql = %d, Nstep_run = %d, Nrun = %d' % (Nstep_eql, Nstep_run, Nrun))


# functions
@jit
def cal_U(rods):
    for i, xi in enumerate(rods):
        for j, xj in enumerate(rods):
            if i != j and np.abs(xi - xj) < 1.0:
                return np.inf
    return 0.0


@jit
def cal_N(rods):
    return rods.size


@jit
def cal_Nc(rods):
    return np.where(np.logical_and(rods >= A, rods <= B))[0].size


@jit
def MC_step(rods):
    rods = shuffle(rods)
    if np.random.rand() < 0.5:
        rods = create(rods)
    else:
        rods = destruct(rods)
    return rods


@jit
def shuffle(rods):
    U0 = cal_U(rods)
    new_rods = np.copy(rods)
    for i, xi in enumerate(rods):
        trial_rods = np.copy(rods)
        dx = (2 * np.random.rand() - 1.0) * dxmax
        # hard wall
        if xi + dx < 0.0:
            trial_rods[i] = -(xi + dx)
        elif xi + dx > Lx:
            trial_rods[i] = 2 * Lx - (xi + dx)
        else:
            trial_rods[i] = xi + dx
        dU = cal_U(trial_rods) - U0
        prob = np.exp(-beta * dU)
        if np.random.rand() < prob:
            new_rods[i] = rods[i]
        rods[i] = xi
    return new_rods


@jit
def create(rods):
    N = cal_N(rods)
    Nc = cal_Nc(rods)
    Vc = B - A
    new_rods = np.copy(rods)
    new_rods.resize(N + 1)
    new_rods[N] = np.random.rand() * Vc + A
    dU = cal_U(new_rods) - cal_U(rods)
    prob = np.exp(-beta * dU + beta * mu) * Vc / (Nc + 1)
    if (np.random.rand() < prob):
        return new_rods
    else:
        return rods


@jit
def destruct(rods):
    Nc = cal_Nc(rods)
    if Nc > 0:
        Vc = B - A
        idx = np.random.choice(np.where(np.logical_and(rods >= A, rods <= B))[0])
        new_rods = np.delete(rods, idx)
        dU = cal_U(new_rods) - cal_U(rods)
        prob = np.exp(-beta * dU - beta * mu) * Nc / Vc
        if np.random.rand() < prob:
            return new_rods
    return rods


def main():
    sumN = 0.0
    sumN2 = 0.0
    sumrho = np.zeros(int(Lx / 0.05))

    for irun in range(Nrun):
        print(irun)
        rods = np.array([Lx / 2.0])

        # equilibrate
        for istep in range(Nstep_eql):
            rods = MC_step(rods)

        # collect data
        for istep in range(Nstep_run):
            rods = MC_step(rods)
            for xi in rods:
                sumrho[int(xi / 0.05)] += 1 / 0.05
            sumN += rods.size
            sumN2 += rods.size**2

    # output
    print('')
    print('# rho dist:')
    for x in sumrho:
        print('%.6f' % (x / Nstep_run / Nrun, ))
    print('')
    print('# <N> = %.6f' % (sumN / Nstep_run / Nrun, ))
    print('# <N2> = %.6f' % (sumN2 / Nstep_run / Nrun, ))



if __name__ == '__main__':
    main()


# END
