import numpy as np
from numba import jit
from matplotlib import pyplot as plt
import seaborn as sns

mass = 1.0
kT = 1.0
beta = 1.0 / kT

dt = 0.001
nu = 0.1
Nstep = 50000
Nstep0 = 20000
N = 500

@jit
def maxwell(N):
    return np.random.normal(0.0, np.sqrt(kT / mass), N)

@jit 
def cal_K(v):
    return np.sum(0.5 * mass * v**2)

@jit
def cal_U(x):
    return np.sum(0.5 * x**2)

@jit
def cal_a(x):
    return -x / mass

@jit
def evolve(x, v):
    a = cal_a(x)
    v += 0.5 * a * dt
    x += v * dt + 0.5 * a * dt**2
    a = cal_a(x)
    v += 0.5 * a * dt
    # Berendsen thermostat
    kTnow = 2 * cal_K(v) / N
    v *= (1.0 + dt / nu * (kT / kTnow - 1.0))**0.5
    '''
    # Andersen thermostat
    randnum = np.random.rand(x.size)
    vnew = maxwell(x.size)
    idx = np.where(randnum < nu*dt)
    v[idx] = vnew[idx]
    '''
    return x, v

def main():
    x = np.random.rand(N)
    v = maxwell(N)
    vrec = np.zeros(N * Nstep)

    for istep in range(1, Nstep0):
        x, v = evolve(x, v)

    U = 0.0
    K = 0.0
    T = 0.0
    for istep in range(1, Nstep):
        U += cal_U(x) / N
        K += cal_K(v) / N
        T += 2 * cal_K(v) / N
        vrec[istep * N:(istep + 1) * N] = v
        print(("%16.6f" * 5) % (istep * dt, U / istep, K / istep, (U+K) / istep, T / istep))
        x, v = evolve(x, v)

    V = np.linspace(-4, 4, 1000)
    P = (beta / 2 / np.pi / mass)**0.5 * np.exp(-beta * V**2 * mass / 2)
    sns.distplot(vrec);
    plt.plot(V, P)
    plt.show()


main()
