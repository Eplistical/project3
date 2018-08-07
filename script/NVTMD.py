import numpy as np
from numba import jit

mass = 1.0
kT = 1.0
dt = 0.01
nu = 5
Nstep = 100000
Nstep0 = 10000
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
    # Andersen thermostat
    randnum = np.random.rand(x.size)
    vnew = maxwell(x.size)
    idx = np.where(randnum < nu*dt)
    v[idx] = vnew[idx]
    return x, v

@jit 
def main():
    x = np.random.rand(N)
    v = maxwell(N)

    U = 0.0
    K = 0.0
    T = 0.0
    for istep in range(1, Nstep):
        U += cal_U(x) / N
        K += cal_K(v) / N
        T += 2 * cal_K(v) / N
        print(("%16.6f" * 5) % (istep * dt, U / istep, K / istep, (U+K) / istep, T / istep))
        x, v = evolve(x, v)

main()
