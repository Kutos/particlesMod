# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pylab as pl
from matplotlib import animation

fig = pl.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('log')
ax.set_xscale('log')

data = np.loadtxt('step_{0}.txt'.format(1))
pl.plot(data[1::,0],data[1::,1],'.-.', label="Initial spectrum")
line, = pl.plot([],[],'.-.', label="Step 1")
pl.xlim(1e-2,1e1)
pl.ylim(1e-5,1e1)
pl.xlabel("kn")
pl.ylabel("Power spectrum")
def init():
    line.set_data([],[])
    return line,
def animate(i):
    data = np.loadtxt('step_{0}.txt'.format(i+1))
    x = data[:,0]
    y = data[:,1]
    line.set_label("Step {0}".format(i+1))
    line.set_data(x[1::],y[1::])
    pl.legend()
    return line,

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=1000, blit=True, interval=20, repeat=False)
ani.save('im.mp4')
pl.show()
