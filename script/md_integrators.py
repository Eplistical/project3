#!/usr/bin/env python3

# compare different MD integrators
# See Frenkel book Ch4 for details

import numpy as np
from numba import jit
import argparse


# config
dt = 0.01
Nstep = int(1e4)
mass = 1.0


def parse_args():
    global dt, Nstep, mass
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nstep', type=int, help="Nstep for running")
    parser.add_argument('--dt', type=float, help="time step")
    parser.add_argument('--mass', type=float, help="mass")
    args = parser.parse_args()
    if args.Nstep is not None:
        Nstep = args.Nstep
    if args.dt is not None:
        dt = args.dt
    if args.mass is not None:
        mass = args.mass

parse_args()
print('# dt = %.4f, Nstep = %d' % (dt, Nstep))
print('# mass = %.4f' % (mass,))

# functions

def cal_U(x):
    """calculate potential energy at position x
    """
    return 0.5 * x**2

def cal_a(x):
    """calculate acceleration according to position x
    """
    return -x / mass


def verlet(x0, v0):
    """Verlet integrator
    """
    xarr = np.zeros(Nstep)
    varr = np.zeros(Nstep)
    xarr[0] = x0
    varr[0] = v0

    # first step
    xlast = x0
    x = x0 + v0 * dt + 0.5 * cal_a(x0) * dt**2

    i = 1
    while i < Nstep:
        xarr[i] = x
        xnext = 2 * x - xlast + cal_a(x) * dt**2
        varr[i] = (xnext - xlast) / 2 / dt

        xlast = x
        x = xnext
        i += 1

    return xarr, varr


def velocity_verlet(x0, v0):
    """velocity Verlet algorithm
    """
    xarr = np.zeros(Nstep)
    varr = np.zeros(Nstep)
    x = x0
    v = v0
    for i in range(0, Nstep):
        xarr[i] = x
        varr[i] = v
        a = cal_a(x)
        x = x + v * dt + 0.5 * a * dt**2
        v = v + 0.5 * (a + cal_a(x)) * dt

    return xarr, varr


def euler(x0, v0):
    xarr = np.zeros(Nstep)
    varr = np.zeros(Nstep)
    x = x0
    v = v0
    for i in range(0, Nstep):
        xarr[i] = x
        varr[i] = v
        a = cal_a(x)
        x = x + v * dt
        v = v + a * dt
    return xarr, varr


def rk4(x0, v0):
    xarr = np.zeros(Nstep)
    varr = np.zeros(Nstep)
    x = x0
    v = v0
    for i in range(0, Nstep):
        xarr[i] = x
        varr[i] = v
        k1, l1 = dt * v,              dt * cal_a(x)
        k2, l2 = dt * (v + 0.5*l1),   dt * cal_a(x + 0.5*k1)
        k3, l3 = dt * (v + 0.5*l2),   dt * cal_a(x + 0.5*k2)
        k4, l4 = dt * (v + l3),       dt * cal_a(x + k3)
        x = x + (k1 + 2*k2 + 2*k3 + k4) / 6
        v = v + (l1 + 2*l2 + 2*l3 + l4) / 6
    return xarr, varr


def Etot(x, v):
    return 0.5*mass*v**2 + cal_U(x)


def main():
    x0 = 2
    v0 = 0
    x1, v1 = euler(x0, v0)
    x2, v2 = verlet(x0, v0)
    x3, v3 = velocity_verlet(x0, v0)
    x4, v4 = rk4(x0, v0)
    for i in range(Nstep):
        print(  i*dt, 
#                x1[i], Etot(x1[i], v1[i]), 
                x2[i], Etot(x2[i], v2[i]), 
                x3[i], Etot(x3[i], v3[i]), 
                x4[i], Etot(x4[i], v4[i])
                )



if __name__ == '__main__':
    main()


# END
