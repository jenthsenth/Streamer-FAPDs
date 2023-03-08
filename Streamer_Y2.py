import numpy as np
import matplotlib.pyplot as plt

def Streamer_Y(T,xs,Re):
    
    # Load in tab delimited file
    filein = 'Streamer_GRL2014_T='+T+'.dat'
    data = np.loadtxt(filein,skiprows=77)
    
    # filein = 'Streamer_WQ_T='+T+'.dat'
    # data = np.loadtxt(filein,skiprows=46)
    # Read in field arrays
    yline = 45
    # yline = 39
        
    yion = data[:,yline]
    nslices = 100
    ndim = yion.size
    
    xmag = data[:,6].reshape((int(ndim/nslices),nslices))
    ymag = - data[:,7].reshape((int(ndim/nslices),nslices))
    vm = data[:,8].reshape((int(ndim/nslices),nslices))
    ftv = data[:,9].reshape((int(ndim/nslices),nslices))
    bmin = data[:,10].reshape((int(ndim/nslices),nslices))
    v = data[:,11].reshape((int(ndim/nslices),nslices))
    birk = data[:,13].reshape((int(ndim/nslices),nslices))
    pedlam = data[:,15].reshape((int(ndim/nslices),nslices))
    pedpsi = data[:,16].reshape((int(ndim/nslices),nslices))
    hall = data[:,17].reshape((int(ndim/nslices),nslices))
    eflux = data[:,18].reshape((int(ndim/nslices),nslices))
    eavg = data[:,19].reshape((int(ndim/nslices),nslices))
    P = data[:,22].reshape((int(ndim/nslices),nslices))
    k = data[:,23].reshape((int(ndim/nslices),nslices))
    Ne = data[:,24].reshape((int(ndim/nslices),nslices))
    Tempe = data[:,25].reshape((int(ndim/nslices),nslices))
    xion = data[:,44].reshape((int(ndim/nslices),nslices))
    yion = - data[:,45].reshape((int(ndim/nslices),nslices))
    
    # xmag = data[:,6].reshape((int(ndim/nslices),nslices))
    # ymag = data[:,7].reshape((int(ndim/nslices),nslices))
    # ftv = data[:,8].reshape((int(ndim/nslices),nslices))
    # bmin = data[:,9].reshape((int(ndim/nslices),nslices))
    # v = data[:,10].reshape((int(ndim/nslices),nslices))
    # birk = data[:,11].reshape((int(ndim/nslices),nslices))
    # P = data[:,12].reshape((int(ndim/nslices),nslices))
    # k = data[:,13].reshape((int(ndim/nslices),nslices))
    # pedlam = data[:,21].reshape((int(ndim/nslices),nslices))
    # hall = data[:,23].reshape((int(ndim/nslices),nslices))
    # eflux = data[:,27].reshape((int(ndim/nslices),nslices))
    # vm = data[:,37].reshape((int(ndim/nslices),nslices))
    # xion = data[:,38].reshape((int(ndim/nslices),nslices))
    # yion = data[:,39].reshape((int(ndim/nslices),nslices))
    
    ny = int(np.sqrt(yion.size))
    nbad = 48
    ny = int(np.sqrt(yion.size)) - nbad
    
    # higher ny causes problems because of chopped values in RCM
    
    xmag = xmag[nbad:]
    ymag = ymag[nbad:]
    vm = vm[nbad:]
    ftv = ftv[nbad:]
    bmin = bmin[nbad:]
    v = v[nbad:]
    birk = birk[nbad:]
    pedlam = pedlam[nbad:]
    pedpsi = pedpsi[nbad:]
    hall = hall[nbad:]
    eflux = eflux[nbad:]
    eavg = eavg[nbad:]
    P = P[nbad:]
    k = k[nbad:]
    Ne = Ne[nbad:]
    Tempe = Tempe[nbad:]
    xion = xion[nbad:]
    yion = yion[nbad:]
    
    c = 3*10**10
    me = 511000/c**2
    e = 1.602*10**(-19)
    dphi = birk*Re
    bi = 5*10**4
    Te = Tempe
    
    eavg_para = Te + e*dphi
    eflux_para = Ne*np.sqrt(Te/(2*np.pi*me))*eavg_para*(bi/bmin - ((bi - bmin)/bmin)*np.exp(-(bmin/(bi - bmin))*e*dphi/Te))
    
    eavg_para = eavg_para/1000
    eflux_para = eflux_para*1.602*10**(-12)*np.where(birk < 0, 1, 0)
    
    Robinson = (40*eavg_para/(16 + eavg_para**2))*np.sqrt(eflux_para)*np.where(birk < 0, 1, 0)
    
    # Plot the RCM-E fields for all slices at the given time step
    # plt.figure(1); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],k[:,i],label='x = ' + str(yion[i,0]))
    # plt.grid(); plt.ylabel('k'); plt.xlabel(r'$x_i$')
    # plt.title(filein)
    # plt.show()
    
    # plt.figure(2); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],ftv[:,i],label='x = ' + str(yion[0,i]))
    # plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$x_i$')
    # plt.title(filein)
    # plt.show()
    
    # plt.figure(3); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],birk[:,i],label='x = ' + str(yion[0,i]))
    # plt.grid(); plt.ylabel(r'$J_\parallel \left[ \mu A/m^2 \right]$'); plt.xlabel(r'$y_i$')
    # plt.title(filein)
    # plt.show()
    
    # plt.figure(4); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],pedpsi[:,i],label='x = ' + str(yion[0,i]))
    # plt.grid(); plt.ylabel(r'$\Sigma_{\phi \phi,diff} \left[ S \right]$'); plt.xlabel(r'$y_i$')
    # plt.title(filein)
    # plt.show()
    
    # plt.figure(5); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],eflux[:,i],label='x = ' + str(yion[0,i]))
    # plt.grid(); plt.ylabel(r'$\bar{\Phi}_\parallel \left[ erg/cm^2 s \right]$'); plt.xlabel(r'$y_i$')
    # plt.title(filein)
    # plt.show()
    
    # plt.figure(6); plt.clf()
    # for i in range(ny):
    #     plt.plot(xion[:,i],Robinson[:,i],label='x = ' + str(yion[0,i]))
    # plt.grid(); plt.ylabel(r'$\Sigma_{\phi \phi,mono} \left[ S \right]$'); plt.xlabel(r'$y_i$')
    # plt.title(filein)
    # plt.show()
    
    # Choose Y-slice
    try:
        iy = int(input('Enter your Y-slice (try ~50, approx. midway): ')) - 1
    except:
        print('Wrong input. Please enter a number ...')
    
    ys = yion[0,iy]
    
    # Define and plot the RCM-E fields for the chosen slice
    x_array = np.flip(xion[:,iy])
    Vm_array = np.flip(vm[:,iy])
    
    Vgrad_array = np.zeros(ny)
    for i in range(ny-1):
        if(i > 0 and i < ny-1):
            Vgrad_array[i] = (Vm_array[i+1] - Vm_array[i-1])/(x_array[i+1] - x_array[i-1])
    
    Vgrad_array[ny - 1] = Vgrad_array[ny - 2]
    Vgrad_array[0] = Vgrad_array[1]
    
    def Vgrad(x):
        y = np.interp(x,x_array,Vgrad_array)
        return y
    
    # Plot alongside analytical profiles and fix parameters
    # Flux tube volume (fixes Scos = S * cos(zeta))
    
    Scos = -3/2*Vgrad(xs)
    Vgradanalytic = -2/3*Scos * np.ones(ny)
    
    plt.figure(7); plt.clf()
    plt.plot(x_array,Vgrad_array,'b-',label = 'RCM-E')
    plt.plot(x_array,Vgradanalytic,'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\nabla V^{-2/3} (x)$'); plt.xlabel(r'$x_i$')
    plt.title(r'FTV Gradient Profile at $y_s$ = ' + str(ys))
    plt.show()
    
    return Scos, ys