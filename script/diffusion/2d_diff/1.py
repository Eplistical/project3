#!/usr/bin/env python3

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns


fig = plt.figure()

#ax = fig.gca()
ax = fig.add_subplot(111, projection='3d')

a = np.loadtxt('1')
#ax.scatter(a[:,1], a[:,2])
shape = [ int(np.sqrt(a.shape[0])) ] * 2
ax.plot_surface(a[:,-3].reshape(shape), a[:,-2].reshape(shape), a[:,-1].reshape(shape))

plt.show()
