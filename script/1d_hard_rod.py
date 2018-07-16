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


def parse_args():
    global A, B, Nstep_eql, Nstep_run
    parser = argparse.ArgumentParser()
    parser.add_argument('--A', type=float, help="control volume begin")
    parser.add_argument('--B', type=float, help="control volume end")
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

parse_args()
print('# Lx = %.2f, A = %.2f, B = %.2f' %(Lx, A, B))
print('# Nstep_eql = %d, Nstep_run = %d' % (Nstep_eql, Nstep_run,))


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
    N = cal_N(rods)
    for i in range(N):
        U0 = cal_U(rods)
        dx = (2 * np.random.rand() - 1.0) * dxmax
        lastxi = rods[i]
        rods[i] += dx
        if rods[i] < 0.0:
            rods[i] *= -1
        elif rods[i] > Lx:
            rods[i] = 2 * Lx - rods[i]
        dU = cal_U(rods) - U0
        prob = np.exp(-beta * dU)
        if np.random.rand() >= prob:
            rods[i] = lastxi
    return rods


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
    print('#' * 32)
    print('# rho dist:')
    for x in sumrho / Nstep_run:
        print(x, end='\t')
    print('')
    print('# <N> = %.6f' % (sumN / Nstep_run, ))
    print('# <N2> = %.6f' % (sumN2 / Nstep_run, ))
    print('#' * 32)



if __name__ == '__main__':
    main()


# END
