#!/usr/bin/env python3

# repeat Sec. II.A results in van Swol's paper 10.1063/1.464945
# 1d lattice w/ GCMC method

import numpy as np
from numba import jit
import argparse

# config
M = 3
Mc = 3
beta = 1.0
mu = -0.45
epsilon = 0.7

Nstep_run = int(10e6)
Nstep_eql = int(5e5)


configurations = [
        [0,0,0],
        [1,0,0],
        [0,1,0],
        [1,1,0],
        [0,0,1],
        [1,0,1],
        [0,1,1],
        [1,1,1],
        ]


# functions
def parse_args():
    global Mc, M
    parser = argparse.ArgumentParser()
    parser.add_argument('--Mc', type=int, help="control volume")
    args = parser.parse_args()
    if args.Mc is not None:
        Mc = args.Mc
    assert Mc <= M


@jit
def cal_N(lattice):
    return np.sum(lattice)


@jit
def cal_Nc(lattice):
    return np.sum(lattice[:Mc])


@jit
def cal_U(lattice):
    return -epsilon * np.sum(lattice[1:] * lattice[:-1])


@jit
def MC_step(lattice):
    shuffle(lattice)
    if np.random.rand() < 0.5:
        create(lattice)
    else:
        destruct(lattice)


@jit
def shuffle(lattice):
    N = cal_N(lattice)
    if N > 0:
        i = np.random.choice(np.where(lattice == 1)[0])
        while True:
            a = np.random.randint(M)
            if i != a:
                break
        if lattice[a] == 0:
            new_lattice = np.copy(lattice)
            new_lattice[i] = 0
            new_lattice[a] = 1
            dU = cal_U(new_lattice) - cal_U(lattice)
            prob = np.exp(-beta * dU)
            if np.random.rand() < prob:
                lattice[i] = 0
                lattice[a] = 1


@jit
def create(lattice):
    idx = np.random.randint(Mc)
    if lattice[idx] == 0:
        new_lattice = np.copy(lattice)
        new_lattice[idx] = 1
        dU = cal_U(new_lattice) - cal_U(lattice)
        Nc = cal_Nc(lattice)
        prob = np.exp(-beta * dU + beta * mu) * Mc / (Nc + 1)
        if np.random.rand() < prob:
            lattice[idx] = 1


@jit
def destruct(lattice):
    Nc = cal_Nc(lattice)
    if Nc > 0:
        idx = np.random.choice(np.where(lattice[:Mc] == 1)[0])
        new_lattice = np.copy(lattice)
        new_lattice[idx] = 0
        dU = cal_U(new_lattice) - cal_U(lattice)
        prob = np.exp(-beta * dU - beta * mu) * Nc / Mc
        if np.random.rand() < prob:
            lattice[idx] = 0


@jit
def main():
    parse_args()
    lattice = np.random.randint(2, size=M)

    # equilibrate
    for istep in range(Nstep_eql):
        MC_step(lattice)

    # collect data
    conf_static = np.zeros(int(2**M))
    for istep in range(Nstep_run):
        MC_step(lattice)
        conf_static[configurations.index(lattice.tolist())] += 1
    conf_static /= Nstep_run

    # post-processing
    rho = np.zeros(M)
    avgN = 0.0
    avgN2 = 0.0
    avg_betaU = 0.0
    for k, conf in zip(conf_static, configurations):
        lattice = np.array(conf)
        rho += k * lattice
        avgN += k * cal_N(lattice)
        avgN2 += k * cal_N(lattice)**2
        avg_betaU += k * -beta * cal_U(lattice)

    # output
    print('# M = %d, Mc = %d' % (M, Mc))
    print('# configuration statistics: ')
    for x in conf_static:
        print(x)
    print('# rho dist:')
    for x in rho:
        print(x)
    print('# <N> = %.6f' % avgN)
    print('# <N2> = %.6f' % avgN2)
    print('# <-beta*U> = %.6f' % avg_betaU)



if __name__ == '__main__':
    main()


# END
