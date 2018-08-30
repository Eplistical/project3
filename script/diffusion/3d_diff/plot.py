#!/usr/bin/env python3

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import sys

fname = sys.argv[1]

Nx = int(sys.argv[2])
layer = int((Nx-1) / 2)

data = np.loadtxt(fname)
# data: [ix, iy, iz, x, y, z, u]

fig, axes = plt.subplots(1, 3, figsize=(18,5))

X = data[:,3].reshape([Nx,Nx,Nx])
Y = data[:,4].reshape([Nx,Nx,Nx])
Z = data[:,5].reshape([Nx,Nx,Nx])
U = data[:,6].reshape([Nx,Nx,Nx])

# z
axes[0].contourf(X[layer,:,:], Y[layer,:,:], U[layer,:,:])
axes[0].set_title('z')
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')

# y
axes[1].contourf(X[:,layer,:], Z[:,layer,:], U[:,layer,:])
axes[1].set_title('y')
axes[1].set_xlabel('x')
axes[1].set_ylabel('z')

# x
axes[2].contourf(Y[:,:,layer], Z[:,:,layer], U[:,:,layer])
axes[2].set_title('x')
axes[2].set_xlabel('y')
axes[2].set_ylabel('z')

plt.show()
