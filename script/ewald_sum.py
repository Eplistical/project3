import numpy as np
from scipy import linalg as la
from scipy.special import erfc
from numba import jit

# implement calculation of Ewald sum


def cal_rhok2(conf, k):
    N = conf.shape[0]
    rhok = 0.0
    for i in range(N):
        ri = conf[i, 0:3]
        qi = conf[i, 3]
        rhok += qi * np.exp(1j * np.sum(k * ri))
    return rhok.real**2 + rhok.imag**2

def cal_ewald_part1(conf, alpha, kimax, L):
    kiarr = np.arange(-kimax, kimax + 1)
    rst = 0.0
    for kx in kiarr:
        for ky in kiarr:
            for kz in kiarr:
                if (kx == 0 and ky == 0 and kz == 0):
                    continue
                else:
                    k = 2 * np.pi / L * np.array([kx, ky, kz])
                    k2 = kx**2 + ky**2 + kz**2
                    rhok2 = cal_rhok2(conf, k)
                    rst += rhok2 / k2 * np.exp(-k2 / alpha / 4)
    return rst * 2 * np.pi / L**3


def cal_ewald_part2(conf, alpha):
    return -np.sqrt(alpha / np.pi) * np.sum(conf[:,3]**2)


def cal_ewald_part3(conf, alpha):
    N = conf.shape[0]
    rst = 0.0
    alpha_sqrt = alpha**0.5
    for i in range(N-1):
        ri = conf[i,0:3]
        qi = conf[i,3]
        for j in range(i+1, N):
            rj = conf[j,0:3]
            qj = conf[j,3]
            rij = la.norm(ri - rj)
            rst += qi * qj * erfc(alpha_sqrt * rij) / rij
    return rst

def cal_ewald_sum(conf, alpha, kimax, L):
    return cal_ewald_part1(conf, alpha, kimax, L) + cal_ewald_part2(conf, alpha) + cal_ewald_part3(conf, alpha)


L = 8.0 # box length
alpha = 1.0 # Ewald parameter
kimax = 8

# charged particle config [x,y,z,q]
conf = np.array([
    [0, 0, 2, -1.0],
    [0, 0, 6, 1.0],
    ])

if __name__ == '__main__':
    print("L = %12.2f" % (L,))
    print("alpha = %12.2f, kimax = %d" % (alpha, kimax))
    print("conf: ")
    print(conf)
    print("Ewald sum: ")
    print(cal_ewald_sum(conf, alpha, kimax, L))
