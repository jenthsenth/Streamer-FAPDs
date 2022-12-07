#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

T = 0 + 60

while T < 61 + 60:

    # Load in tab delimited file
    filein = 'Streamer_WQ_T='+str(T)+'.dat'
    data = np.loadtxt(filein,skiprows=77)
    
    # Read in field arrays
    # xline = 44
    xline = 38
    blah = data[:,xline]
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
    
    # Plot the RCM-E fields for all slices at the given time step
    plt.figure(1); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],pedlam[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel(r'$\Sigma(y_i)$'); plt.xlabel(r'$y_i$')
    plt.title('Ionospheric Pedersen Conductivity at T = '+str(T))
    plt.show()
    
    T += 1

# T = 0

# while T < 31:

#     # Load in tab delimited file
#     filein = 'Streamer_GRL2014_T='+str(T)+'.dat'
#     data = np.loadtxt(filein,skiprows=77)
    
#     # Read in field arrays
#     yion = data[:,45]
#     nslices = 100
#     ndim = yion.size
    
#     xmag = data[:,6].reshape((int(ndim/nslices),nslices))
#     ymag = data[:,7].reshape((int(ndim/nslices),nslices))
#     vm = data[:,8].reshape((int(ndim/nslices),nslices))
#     ftv = data[:,9].reshape((int(ndim/nslices),nslices))
#     bmin = data[:,10].reshape((int(ndim/nslices),nslices))
#     v = data[:,11].reshape((int(ndim/nslices),nslices))
#     birk = data[:,13].reshape((int(ndim/nslices),nslices))
#     pedlam = data[:,15].reshape((int(ndim/nslices),nslices))
#     hall = data[:,17].reshape((int(ndim/nslices),nslices))
#     eflux = data[:,18].reshape((int(ndim/nslices),nslices))
#     N = data[:,20].reshape((int(ndim/nslices),nslices))
#     Temp = data[:,21].reshape((int(ndim/nslices),nslices))
#     P = data[:,22].reshape((int(ndim/nslices),nslices))
#     k = data[:,23].reshape((int(ndim/nslices),nslices))
#     Ne = data[:,24].reshape((int(ndim/nslices),nslices))
#     Tempe = data[:,25].reshape((int(ndim/nslices),nslices))
#     Pe = data[:,26].reshape((int(ndim/nslices),nslices))
#     ke = data[:,27].reshape((int(ndim/nslices),nslices))
#     Ni = data[:,28].reshape((int(ndim/nslices),nslices))
#     Tempi = data[:,20].reshape((int(ndim/nslices),nslices))
#     Pi = data[:,30].reshape((int(ndim/nslices),nslices))
#     ki = data[:,31].reshape((int(ndim/nslices),nslices))
#     xion = data[:,44].reshape((int(ndim/nslices),nslices))
#     yion = data[:,45].reshape((int(ndim/nslices),nslices))
#     brat = data[:,51].reshape((int(ndim/nslices),nslices))
#     edistrat = data[:,52].reshape((int(ndim/nslices),nslices))
#     btotrat = data[:,57].reshape((int(ndim/nslices),nslices))
    
#     ny = int(np.sqrt(yion.size))
    
#     # Plot the RCM-E fields for all slices at the given time step
#     plt.figure(1); plt.clf()
#     for i in range(ny):
#         plt.plot(xion[:,i],edistrat[:,i],label='x = ' + str(yion[i,0]))
#     plt.grid(); plt.ylabel('edistrat'); plt.xlabel(r'$x_i$')
#     plt.title(filein)
#     plt.show()
    
#     T += 1