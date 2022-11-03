import numpy as np
import matplotlib.pyplot as plt

def Streamer_Y(T,xs):
    
    # T = str(10)
    # xs = -0.3969697058
    # Load in tab delimited file
    filein = 'Streamer_GRL2014_T='+T+'.dat'
    data = np.loadtxt(filein,skiprows=77)
    
    # Read in field arrays
    yion = data[:,45]
    nslices = 100
    ndim = yion.size
    
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
    
    ny = int(np.sqrt(yion.size))
    nbadL = 48
    ny = int(np.sqrt(yion.size)) - nbadL
    
    # higher ny causes problems because of chopped values in RCM
    
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
    for i in range(ny):
        plt.plot(xion[:,i],k[:,i],label='x = ' + str(yion[i,0]))
    plt.grid(); plt.ylabel('k'); plt.xlabel(r'$x_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(2); plt.clf()
    for i in range(ny):
        plt.plot(xion[:,i],ftv[:,i],label='x = ' + str(yion[0,i]))
    plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$x_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(3); plt.clf()
    for i in range(ny):
        plt.plot(xion[:,i],pedlam[:,i],label='x = ' + str(yion[i,0]))
    plt.grid(); plt.ylabel('pedlam'); plt.xlabel(r'$x_i$')
    plt.title(filein)
    plt.show()
    
    # Choose Y-slice
    try:
        iy = int(input('Enter your Y-slice: ')) - 1
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
    
    def Vgrad(x):
        y = np.interp(x,x_array,Vgrad_array)
        return y
    
    plt.figure(7); plt.clf()
    plt.plot(x_array,Vgrad_array)
    plt.grid(); plt.ylabel(r'$\nabla V^{-2/3} (x)$'); plt.xlabel(r'$x_i$')
    plt.title(r'Vgrad Term at $T = $' + str(T) + ' and $Y = $' + str(ys))
    plt.show()
    
    # Plot alongside analytical profiles and fix parameters
    # Flux tube volume (fixes Scos = S * cos(zeta))
    
    Scos = -3/2*Vgrad(xs)
    Vgradanalytic = -2/3*Scos * np.ones(ny)
    
    plt.figure(8); plt.clf()
    plt.plot(x_array,Vgrad_array,'b-',label = 'RCM-E')
    plt.plot(x_array,Vgradanalytic,'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\nabla V^{-2/3} (x)$'); plt.xlabel(r'$x_i$')
    plt.title(r'FTV Gradient Profile at $y_s$ = ' + str(ys))
    plt.show()
    
    return Scos, ys