#!/usr/bin/env python3

# potentials for H & H+

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt


def LJ(r, sigma, epsilon):
    ir2 = (sigma / r)**2
    ir6 = ir2**3
    U = 4.0 * epsilon * ir6 * (ir6 - 1.0)
    W = 16.0 * epsilon * ir6 * (ir6 - 0.5)
    return U



sigma = 5.0
epsilon = 2.0
F = 0.0
z = np.linspace(sigma*0.8, 20, 1000)
plt.plot(z, LJ(z, sigma, epsilon) + F, label='Solvant')

sigma = 3.0
epsilon = 4.0
F = 0.0
z = np.linspace(sigma*0.8, 20, 1000)
plt.plot(z, LJ(z, sigma, epsilon) + F, label='H')

sigma = 4.0
epsilon = 5.0
F = 0.0
z = np.linspace(sigma*0.8, 20, 1000)
plt.plot(z, LJ(z, sigma, epsilon) + F, label='H+')

plt.ylim([-10, 50])
plt.legend()
plt.show()
