#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 10:16:39 2022

@author: toffo
"""
import numpy as np
import matplotlib.pyplot as plt

file='Streamer File.dat'
data = np.loadtxt(file,skiprows=65)

x1 = data[:,44]
y1 = data[:,45]
pvg = data[:,23]

nslices = 100
ndim = x1.size
x1=x1.reshape((int(ndim/nslices),nslices)) # row const
y1=y1.reshape((int(ndim/nslices),nslices)) # column const
pvg=pvg.reshape((int(ndim/nslices),nslices))

# data.reshape((n/ncols,ncols))

xindex = 50


plt.figure(1)
plt.clf()


plt.plot(y1[xindex,:],pvg[xindex,:])
plt.title(' for x ='+str(x1[xindex,0]))


plt.xlabel('y')
plt.ylabel(r'$pV^\gamma$')
# plt.draw()
plt.show()