# -*- coding: utf-8 -*-
import numpy as np
import pylab as pl

fig = pl.figure(figsize=(10, 10))
ax = fig.add_subplot(2, 1, 1)
ax.set_yscale('log')
ax.set_xscale('log')

data = np.loadtxt('pk_step_10.txt')
pl.plot(data[:,0],data[:,1], '.', label="Final power spectrum")

data = np.loadtxt('initial_pk_linear_theo.txt')
ax.plot(data[:,0],data[:,1], '.', label="Initial power spectrum")

pl.title("Power Specturm")
pl.xlabel("kn")
pl.ylabel("Pk")

pl.legend()
pl.savefig('pk_step_1000.pdf', dpi=fig.dpi)