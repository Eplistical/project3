import sys
import numpy as np
from scipy import integrate

Nx = 40
dx = 0.05
D = 1.0
dt = 0.0001

assert 2 * D * dt / dx**2 < 1.0

Nstep = int(sys.argv[1])
Dinvdx2= D / dx**2
Ntot = Nx * Nx

def cal_dudt(t,u):
    dudt = np.zeros(Ntot)
    for ix in range(1,Nx-1):
        for iy in range(1,Nx-1):
            idx = ix + Nx * iy
            dudt[idx] = Dinvdx2 * (u[idx+1] + u[idx-1] + u[idx+Nx] + u[idx-Nx] - 4*u[idx])
    return dudt


def make_jac(t, u):
    jac = np.zeros([Ntot, Ntot])
    for ix in range(1,Nx-1):
        for iy in range(1,Nx-1):
            idx = ix + Nx * iy
            jac[idx, idx] = -4 * Dinvdx2 
            jac[idx, idx+1] = Dinvdx2 
            jac[idx, idx-1] = Dinvdx2 
            jac[idx, idx+Nx] = Dinvdx2 
            jac[idx, idx-Nx] = Dinvdx2 
    return jac


def init_u():
    mu = Nx * dx / 2
    sigma = 0.05
    print(0.5 / np.pi / sigma**2)

    u = np.zeros(Ntot)
    for ix in range(1,Nx-1):
        x = ix * dx
        for iy in range(1,Nx-1):
            y = iy * dx
            idx = ix + Nx * iy
            u[idx] = 0.5 / np.pi / sigma**2 * np.exp(-0.5 * (((x-mu) / sigma)**2 + ((y-mu) / sigma)**2))
    return u

rr = integrate.ode(cal_dudt, make_jac).set_integrator('vode', method='bdf', atol='1e-14', rtol='1e-4')
rr.set_initial_value(init_u(), 0.0)
Anastep = 10

for i in range(Nstep):
    print("step ", i, "t = ", rr.t)
    if i % Anastep == 0:
        with open('p%d' % (i / Anastep), 'w') as f:
            for ix in range(Nx):
                for iy in range(Nx):
                    f.write('%16.6f%16.6f%16.6f\n' % (ix*dx, iy*dx, rr.y[ix + Nx * iy]))

    rr.integrate(rr.t + dt)

    # boundary
    for iy in range(Nx):
        rr.y[0 + iy * Nx] = 0.0
        rr.y[Nx-1 + iy * Nx] = 0.0
    for ix in range(Nx):
        rr.y[ix + 0 * Nx] = 0.0
        rr.y[ix + (Nx-1) * Nx] = 0.0

    rr.set_initial_value(rr.y, rr.t)

