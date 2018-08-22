#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

D = 1.0
dt = 0.001
Nx = 501
Nt = 2000
dt = 0.0001


x, dx = np.linspace(-4, 4, Nx, retstep=True)
print("(accurate when < 1) 2Ddt/dx**2 = ", 2*D*dt/dx**2)
assert(2*D*dt/dx**2 < 1.0)

def init_u(x):
    mu = 0.0
    sigma = 0.05
    u = np.exp(-(x-mu)**2 / 2 / sigma**2) / np.sqrt(2 * np.pi * sigma**2)
    u[1] = 0.0
    u[-1] = 0.0
    return u


def cal_dudt(u, dx, D):
    """calc du/dt = D * d^2u/dx^2"""
    N = u.size
    dudt = np.zeros(N)
    for i in range(1,N-1):
        dudt[i] = (u[i+1] + u[i-1] - 2 * u[i]) / dx**2
    dudt[0] = 0.0
    dudt[-1] = 0.0
    return D * dudt

fig, ax = plt.subplots()

u = init_u(x)

for istep in range(Nt):
    # RK4
    k1 = dt * cal_dudt(u, dx, D)
    k2 = dt * cal_dudt(u + 0.5 * k1, dx, D)
    k3 = dt * cal_dudt(u + 0.5 * k2, dx, D)
    k4 = dt * cal_dudt(u + k3, dx, D)
    u = u + (k1 + 2*k2 + 2*k3 + k4) / 6.0 
    if istep % int(Nt / 10) == 0:
        ax.plot(x, u, label="%f" % (istep * dt, ))

plt.legend()
plt.show()
