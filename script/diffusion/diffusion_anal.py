#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

def cal_u(x, t, D = 1.0):
    assert t > 0.0
    return (4*np.pi*D*t)**-0.5 * np.exp(-x**2/4/D/t)


D = 1.0
N = 1000
x = np.linspace(-10, 10, N)

fig, ax = plt.subplots()


for i, t in enumerate(np.arange(1, 50) * 0.05):
    print(i+1)
    ax.cla()
    ax.plot(x, cal_u(x, t, D))
    ax.set_xlim([-10,10])
    ax.set_ylim([0,1.5])
    fig.savefig('a%d.jpg' % (i+1))
