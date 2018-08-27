import numpy as np
from scipy import integrate

Nx = 100
Nstep = 500
D = 1.0
G0 = 0.0
V0 = 1.0
V1 = -1.0
scan_rate = 1.0
dt = abs(V0 - V1) / (Nstep * scan_rate)
beta = 38.92
alpha = 0.2
nFA = 1.0
kA = 0.0
k0 = 1.0

kf = k0 * np.exp(-beta * alpha * (V0 - G0))
kb = k0 * np.exp(beta * (1 - alpha) * (V0 - G0))

dx = 6 * np.sqrt(2. * abs(V0 - V1) * D / scan_rate) / Nx

Dinvdx2 = D / dx / dx
delta = 1 / dx
cAinf = 1.0
cBinf = 0.0


def cal_dydt(t,y):
    dydt = np.zeros(Nx*2)

    dydt[0] = Dinvdx2 * (y[1] - y[0]) + kA * y[Nx] - (kf * y[0] - kb * y[Nx]) * delta
    dydt[Nx] = Dinvdx2 * (y[1+Nx] - y[Nx])- kA * y[Nx] + (kf * y[0] - kb * y[Nx]) * delta

    for i in range(1,Nx-1):
        dydt[i] = Dinvdx2 * (y[i-1] - 2 * y[i] + y[i+1]) + kA * y[i+Nx] 
        dydt[i+Nx] = Dinvdx2 * (y[i-1+Nx] - 2 * y[i+Nx] + y[i+1+Nx]) - kA * y[i+Nx] 

    dydt[Nx-1] = Dinvdx2 * (y[Nx-1] - 2 * y[Nx-2] + y[Nx-3]) + kA * y[Nx-1+Nx]
    dydt[Nx-1+Nx] = Dinvdx2 * (y[Nx-1+Nx] - 2 * y[Nx-2+Nx] + y[Nx-3+Nx]) - kA * y[Nx-1] 

    return dydt


def make_jac(t, y):
    jac = np.zeros((Nx*2, Nx*2))

    jac[0,0] = -Dinvdx2 - kf * delta
    jac[0,1] = Dinvdx2
    jac[0,Nx] = kA + kb * delta

    jac[Nx-1,Nx-3] = Dinvdx2
    jac[Nx-1,Nx-2] = -2 * Dinvdx2
    jac[Nx-1,Nx-1] = Dinvdx2
    jac[Nx-1,Nx-1] = kA 

    jac[Nx,Nx] = -Dinvdx2 - kb * delta - kA
    jac[Nx,Nx+1] = Dinvdx2
    jac[Nx,0] = kf * delta

    jac[Nx-1+Nx,Nx-3+Nx] = Dinvdx2
    jac[Nx-1+Nx,Nx-2+Nx] = -2 * Dinvdx2
    jac[Nx-1+Nx,Nx-1+Nx] = Dinvdx2

    for i in range(1,Nx-1):
        jac[i,i] = -2 * Dinvdx2
        jac[i,i+Nx] = kA
        jac[i+Nx,i+Nx] = -2 * Dinvdx2 - kA

        jac[i,i+1] = Dinvdx2
        jac[i,i-1] = Dinvdx2

        jac[i+Nx,i+1+Nx] = Dinvdx2
        jac[i+Nx,i-1+Nx] = Dinvdx2

    return jac


rr = integrate.ode(cal_dydt, make_jac).set_integrator('vode', method='bdf', atol='1e-14', rtol='1e-4')
rr.set_initial_value([cAinf] * Nx + [cBinf] * Nx, 0.0)

V = V0
for i in range(Nstep):
    V -= dt

    kf = k0 * np.exp(-beta * alpha * (V - G0))
    kb = k0 * np.exp(beta * (1 - alpha) * (V - G0))

    rr.integrate(rr.t + dt)

    # x = inf boundary
    rr.y[Nx-1] = cAinf
    rr.y[2*Nx-1] = cBinf

    rr.set_initial_value(rr.y, rr.t)

    # measure current
    print('%20.10f%20.10f' % (V, nFA * (kf * rr.y[0] - kb * rr.y[Nx])))

# Return loop
V = V1
for i in range(Nstep):
    V += dt

    kf = k0 * np.exp(-beta * alpha * (V - G0))
    kb = k0 * np.exp(beta * (1 - alpha) * (V - G0))

    rr.integrate(rr.t + dt)

    # x = inf boundary
    rr.y[Nx-1] = cAinf
    rr.y[2*Nx-1] = cBinf

    rr.set_initial_value(rr.y, rr.t)

    # measure current
    print('%20.10f%20.10f' % (V, nFA * (kf * rr.y[0] - kb * rr.y[Nx])))
