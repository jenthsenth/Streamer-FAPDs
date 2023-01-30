#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

T = 0

while T < 31:

    # Load in tab delimited file
    filein = 'Streamer_GRL2014_T='+str(T)+'.dat'
    data = np.loadtxt(filein,skiprows=77)
    
    # Read in field arrays
    blah = data[:,44]
    nslices = 100
    ndim = blah.size
    
    xmag = data[:,6].reshape((int(ndim/nslices),nslices))
    ymag = data[:,7].reshape((int(ndim/nslices),nslices))
    vm = data[:,8].reshape((int(ndim/nslices),nslices))
    ftv = data[:,9].reshape((int(ndim/nslices),nslices))
    bmin = data[:,10].reshape((int(ndim/nslices),nslices))
    v = data[:,11].reshape((int(ndim/nslices),nslices))
    birk = data[:,13].reshape((int(ndim/nslices),nslices))
    pedlam = data[:,15].reshape((int(ndim/nslices),nslices))
    hall = data[:,17].reshape((int(ndim/nslices),nslices))
    eflux = data[:,18].reshape((int(ndim/nslices),nslices))
    N = data[:,20].reshape((int(ndim/nslices),nslices))
    Temp = data[:,21].reshape((int(ndim/nslices),nslices))
    P = data[:,22].reshape((int(ndim/nslices),nslices))
    k = data[:,23].reshape((int(ndim/nslices),nslices))
    Ne = data[:,24].reshape((int(ndim/nslices),nslices))
    Tempe = data[:,25].reshape((int(ndim/nslices),nslices))
    Pe = data[:,26].reshape((int(ndim/nslices),nslices))
    ke = data[:,27].reshape((int(ndim/nslices),nslices))
    Ni = data[:,28].reshape((int(ndim/nslices),nslices))
    Tempi = data[:,20].reshape((int(ndim/nslices),nslices))
    Pi = data[:,30].reshape((int(ndim/nslices),nslices))
    ki = data[:,31].reshape((int(ndim/nslices),nslices))
    xion = data[:,44].reshape((int(ndim/nslices),nslices))
    yion = data[:,45].reshape((int(ndim/nslices),nslices))
    brat = data[:,51].reshape((int(ndim/nslices),nslices))
    edistrat = data[:,52].reshape((int(ndim/nslices),nslices))
    btotrat = data[:,57].reshape((int(ndim/nslices),nslices))
    
    # Excise points outside of the boundary
    nbadL = 61
    nx = int(np.sqrt(xion.size)) - nbadL
    
    xmag = xmag[nbadL:]
    ymag = ymag[nbadL:]
    vm = vm[nbadL:]
    ftv = ftv[nbadL:]
    bmin = bmin[nbadL:]
    v = v[nbadL:]
    birk = birk[nbadL:]
    pedlam = pedlam[nbadL:]
    hall = hall[nbadL:]
    eflux = eflux[nbadL:]
    N = N[nbadL:]
    Temp = Temp[nbadL:]
    P = P[nbadL:]
    k = k[nbadL:]
    Ne = Ne[nbadL:]
    Tempe = Tempe[nbadL:]
    Pe = Pe[nbadL:]
    ke = ke[nbadL:]
    Ni = Ni[nbadL:]
    Tempi = Tempi[nbadL:]
    Pi = Pi[nbadL:]
    ki = ki[nbadL:]
    xion = xion[nbadL:]
    yion = yion[nbadL:]
    brat = brat[nbadL:]
    edistrat = edistrat[nbadL:]
    btotrat = btotrat[nbadL:]
    
    xion = xion[:,0]
    yion = yion[0,:]
    
    # Plot the RCM-E fields for all slices at the given time step
    
    field = eflux
    X, Y = np.meshgrid(yion, xion)
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, field)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Eflux at T = '+str(T))
    ax.set_ylabel(r'$x_i$')
    ax.set_xlabel(r'$y_i$')
    
    # calc index of min/max Z value
    # xmin, ymin = np.unravel_index(np.argmin(field), field.shape)
    # xmax, ymax = np.unravel_index(np.argmax(field), field.shape)
    
    # ax.plot(yion[ymin], xion[xmin], 'bo')
    # ax.plot(yion[ymax], xion[xmax], 'ro')
    
    plt.show()
    
    T += 1