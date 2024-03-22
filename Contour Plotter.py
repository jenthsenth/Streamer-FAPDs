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
    eavg = data[:,19].reshape((int(ndim/nslices),nslices))
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
    nbad = 61
    nx = int(np.sqrt(xion.size)) - nbad
    
    xmag = xmag[nbad:]
    ymag = ymag[nbad:]
    vm = vm[nbad:]
    ftv = ftv[nbad:]
    bmin = bmin[nbad:]
    v = v[nbad:]
    birk = birk[nbad:]
    pedlam = pedlam[nbad:]
    hall = hall[nbad:]
    eflux = eflux[nbad:]
    eavg = eavg[nbad:]
    N = N[nbad:]
    Temp = Temp[nbad:]
    P = P[nbad:]
    k = k[nbad:]
    Ne = Ne[nbad:]
    Tempe = Tempe[nbad:]
    Pe = Pe[nbad:]
    ke = ke[nbad:]
    Ni = Ni[nbad:]
    Tempi = Tempi[nbad:]
    Pi = Pi[nbad:]
    ki = ki[nbad:]
    xion = xion[nbad:]
    yion = yion[nbad:]
    brat = brat[nbad:]
    edistrat = edistrat[nbad:]
    btotrat = btotrat[nbad:]
    
    xion = xion[:,0]
    yion = yion[0,:]
    
    c = 3*10**10
    me = 511000/c**2
    e = 1.602*10**(-19)
    dphi = birk*2500
    bi = 5*10**4
    Te = Tempe
    
    eavg_para = Te + e*dphi
    eflux_para = Ne*np.sqrt(Te/(2*np.pi*me))*eavg_para*(bi/bmin - ((bi - bmin)/bmin)*np.exp(-(bmin/(bi - bmin))*e*dphi/Te))
    
    eavg_para = eavg_para/1000
    eflux_para = eflux_para*1.602*10**(-12)*np.where(birk < 0, 1, 0)
    
    Robinson = (40*eavg_para/(16 + eavg_para**2))*np.sqrt(eflux_para)
    
    # Plot the RCM-E fields for all slices at the given time step
    
    field = k
    X, Y = np.meshgrid(yion, xion)
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, field)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Ionospheric Entropy at T = '+str(T))
    ax.set_ylabel(r'$x_i$')
    ax.set_xlabel(r'$y_i$')
    
    # # calc index of min/max Z value
    xmin, ymin = np.unravel_index(np.argmin(field), field.shape)
    xmax, ymax = np.unravel_index(np.argmax(field), field.shape)
    
    ax.plot(yion[ymin], xion[xmin], 'bo')
    ax.plot(yion[ymax], xion[xmax], 'ro')
    
    plt.show()
    
    T += 1