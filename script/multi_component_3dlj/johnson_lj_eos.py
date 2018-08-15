#!/usr/bin/env python3
import numpy as np
import sys
from matplotlib import pyplot as plt

# --- 32 linear para + 1 nonlinear para --- #
x1 =  0.8623085097507421
x2 =  2.976218765822098
x3 = -8.402230115796038
x4 =  0.1054136629203555
x5 = -0.8564583828174598
x6 =  1.582759470107601
x7 =  0.7639421948305453
x8 =  1.753173414312048
x9 =  2.798291772190376E+03
x10 = -4.8394220260857657E-02
x11 =  0.9963265197721935
x12 = -3.698000291272493E+01
x13 =  2.084012299434647E+01
x14 =  8.305402124717285E+01
x15 = -9.574799715203068e+02
x16 = -1.477746229234994E+02
x17 =  6.398607852471505E+01
x18 =  1.603993673294834E+01
x19 =  6.805916615864377E+01
x20 = -2.791293578795945E+03
x21 = -6.245128304568454
x22 = -8.116836104958410E+03
x23 =  1.488735559561229E+01
x24 = -1.059346754655084E+04
x25 = -1.131607632802822E+02
x26 = -8.867771540418822E+03
x27 = -3.986982844450543E+01
x28 = -4.689270299917261E+03
x29 =  2.593535277438717E+02
x30 = -2.694523589434903E+03
x31 = -7.218487631550215E+02
x32 =  1.721802063863269E+02
gm = 3.0
# --- 32 linear para + 1 nonlinear para --- #

def cal_a(T):
    a1 = x1*T +  x2*T**0.5 +  x3 +  x4*T**-1 +  x5*T**-2
    a2 = x6*T              +  x7 +  x8*T**-1 +  x9*T**-2
    a3 = x10*T             +  x11 + x12*T**-1
    a4 =                      x13
    a5 =                            x14*T**-1 + x15*T**-2
    a6 =                            x16*T**-1
    a7 =                            x17*T**-1 + x18*T**-2
    a8 =                                        x19*T**-2
    return a1,a2,a3,a4,a5,a6,a7,a8


def cal_b(T):
    b1 = x20*T**-2 +  x21*T**-3
    b2 = x22*T**-2              +  x23*T**-4
    b3 = x24*T**-2 +  x25*T**-3
    b4 = x26*T**-2              +  x27*T**-4
    b5 = x28*T**-2 +  x29*T**-3
    b6 = x30*T**-2 +  x31*T**-3 +  x32*T**-4
    return b1,b2,b3,b4,b5,b6


def cal_c(T):
    c1 = 0.5*x2*T**0.5 +  x3 +  2*x4*T**-1 +  3*x5*T**-2
    c2 =                  x7 +  2*x8*T**-1 +  3*x9*T**-2
    c3 =                  x11 + 2*x12*T**-1 
    c4 =                  x13
    c5 =                        2*x14*T**-1 + 3*x15*T**-2
    c6 =                        2*x16*T**-1
    c7 =                        2*x17*T**-1 + 3*x18*T**-2
    c8 =                                      3*x19*T**-2
    return c1,c2,c3,c4,c5,c6,c7,c8


def cal_d(T):
    d1 = 3*x20*T**-2 +  4*x21*T**-3
    d2 = 3*x22*T**-2                +  5*x23*T**-4
    d3 = 3*x24*T**-2 +  4*x25*T**-3
    d4 = 3*x26*T**-2                +  5*x27*T**-4
    d5 = 3*x28*T**-2 +  4*x29*T**-3
    d6 = 3*x30*T**-2 +  4*x31*T**-3 +  5*x32*T**-4
    return d1,d2,d3,d4,d5,d6


def cal_G(rho):
    F = np.exp(-gm*rho*rho)
    G1 = (1 - F) / 2 / gm
    G2 = -(F*rho**2 - 2*G1) / 2 / gm
    G3 = -(F*rho**4 - 4*G2) / 2 / gm
    G4 = -(F*rho**6 - 6*G3) / 2 / gm
    G5 = -(F*rho**8 - 8*G4) / 2 / gm
    G6 = -(F*rho**10 - 10*G5) / 2 / gm
    return G1,G2,G3,G4,G5,G6


def cal_P(rho, T):
    F = np.exp(-gm*rho*rho)
    a1,a2,a3,a4,a5,a6,a7,a8 = cal_a(T)
    b1,b2,b3,b4,b5,b6 = cal_b(T)
    P = rho * T
    P += a1*rho**2 + a2*rho**3 + a3*rho**4 + a4*rho**5 \
        + a5*rho**6 + a6*rho**7 + a7*rho**8 + a8*rho**9
    P += F * (b1*rho**3 + b2*rho**5 + b3*rho**7 \
        + b4*rho**9 + b5*rho**11 + b6*rho**13)
    return P

def cal_U(rho, T):
    c1,c2,c3,c4,c5,c6,c7,c8 = cal_c(T)
    d1,d2,d3,d4,d5,d6 = cal_d(T)
    G1,G2,G3,G4,G5,G6 = cal_G(rho)
    U = 0.0
    U += (c1*rho**1)/1 + (c2*rho**2)/2 + (c3*rho**3)/3 + (c4*rho**4)/4  \
            + (c5*rho**5)/5 + (c6*rho**6)/6 + (c7*rho**7)/7 + (c8*rho**8)/8
    U += d1*G1 + d2*G2 + d3*G3 + d4*G4 + d5*G5 + d6*G6
    return U

def cal_A(rho, T):
    a1,a2,a3,a4,a5,a6,a7,a8 = cal_a(T)
    b1,b2,b3,b4,b5,b6 = cal_b(T)
    G1,G2,G3,G4,G5,G6 = cal_G(rho)

    A = 0.0
    A += (a1*rho**1)/1 + (a2*rho**2)/2 + (a3*rho**3)/3 + (a4*rho**4)/4  \
            + (a5*rho**5)/5 + (a6*rho**6)/6 + (a7*rho**7)/7 + (a8*rho**8)/8
    A += b1*G1 + b2*G2 + b3*G3 + b4*G4 + b5*G5 + b6*G6

    return A

def cal_mu(rho, T):
    A = cal_A(rho, T)
    P = cal_P(rho, T)
    mu_ex = A + P / rho - T
    mu_0 = T * np.log(rho)
    mu = mu_0 + mu_ex
    return mu

if len(sys.argv) > 2:
    try:
        print("%10s%16s%16s%16s%16s" % ("T", "rho", "U", "P", "mu"))
        i = 0
        while True:
            T = float(sys.argv[2*i+1])
            rho = float(sys.argv[2*i+2])
            P = cal_P(rho,T)
            U = cal_U(rho,T)
            mu = cal_mu(rho, T)
            print("%10.2f%16.8f%16.8f%16.8f%16.8f" % (T, rho, U, P, mu))
            i += 1
    except:
        sys.exit()


data = np.loadtxt('johnson_data.txt')
Tarr = data[:,0]
rhovaparr = data[:,1]
rholiqarr = data[:,2]
Parr = data[:,3]

Pvap = np.zeros(Tarr.size)
Pliq = np.zeros(Tarr.size)

for i, T, rhovap, rholiq, P in zip(np.arange(Tarr.size), Tarr, rhovaparr, rholiqarr, Parr):
    Pvap[i] = cal_P(rhovap, T)
    Pliq[i] = cal_P(rholiq, T)
    print("%10.2f%16.8f%16.8f" % (T, rhovap, cal_U(rhovap, T)))

"""
plt.plot(rhovaparr, Pvap, ':b')
plt.plot(rholiqarr, Pliq, ':b')
plt.plot(rhovaparr, Parr, '-r')
plt.plot(rholiqarr, Parr, '-r')
"""

plt.plot(Tarr, Pvap, ':b')
plt.plot(Tarr, Pliq, ':g')
plt.plot(Tarr, Parr, '-r')
plt.show()

# END
