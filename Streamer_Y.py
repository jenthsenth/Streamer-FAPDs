import numpy as np
import matplotlib.pyplot as plt

def Streamer_Y(T,xspecial):
    
    # Load in tab delimited file
    filein = 'T'+T+'yi.txt'
    data = np.loadtxt(filein,skiprows=1)
    
    # Read in field arrays
    xin = data[:,38]
    yin = data[:,39]
    ftv = data[:,8]
    bmin = data[:,9]
    v = data[:,10]
    birk = data[:,11]
    p = data[:,12]
    k = data[:,13]
    vtotx = data[:,14]
    vtoty = data[:,15]
    ey = data[:,20]
    pedlam = data[:,21]
    hall = data[:,23]
    tp = data[:,24]
    te = data[:,25]
    ne = data[:,26]
    eflux = data[:,27]
    eavg = data[:,28]
    vm = data[:,37]
    
    # Find unique values of y
    yopts = yin[0]
    ylast = yopts
    
    # Count unique Y-slices and determine corresponding x-indices
    ix = np.array(0)
    
    for i in range(yin.size):
        if(yin[i]!=ylast):
            ix = np.append(ix,i)
            yopts = np.append(yopts,yin[i])
            ylast = yin[i]
    
    ix = np.append(ix,yin.size)
    iy = yopts.size
    
    print('There are ',iy,' unique y values in this dataset')
    print('y = ',yopts)
    
    # Plot the RCM-E fields for all slices at the given time step
    ixmin = ix[0:iy]
    ixmax = ix[1:]
    
    plt.figure(1); plt.clf()
    for i in range(yopts.size):
        plt.plot(xin[ixmin[i]:ixmax[i]],k[ixmin[i]:ixmax[i]],label='y = ' + str(yopts[i]))
    plt.legend(); plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$x$')
    plt.title(filein)
    plt.show()
    
    plt.figure(2); plt.clf()
    for i in range(yopts.size):
        plt.plot(xin[ixmin[i]:ixmax[i]],ftv[ixmin[i]:ixmax[i]],label='y = ' + str(yopts[i]))
    plt.legend(); plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$x$')
    plt.title(filein)
    plt.show()
    
    # plt.figure(3); plt.clf()
    # for i in range(yopts.size):
    #     plt.plot(xin[ixmin[i]:ixmax[i]],vm[ixmin[i]:ixmax[i]],label='y = ' + str(yopts[i]))
    # plt.legend(); plt.grid(); plt.ylabel('vm'); plt.xlabel(r'$x$')
    # plt.title(filein)
    # plt.show()
    
    plt.figure(3); plt.clf()
    for i in range(yopts.size):
        plt.plot(xin[ixmin[i]:ixmax[i]],pedlam[ixmin[i]:ixmax[i]],label='y = ' + str(yopts[i]))
    plt.legend(); plt.grid(); plt.ylabel('pedlam'); plt.xlabel(r'$x$')
    plt.title(filein)
    plt.show()
    
    # Choose Y-slice
    try:
        iy = int(input('Enter your Y-slice: ')) - 1
    except:
        print('Wrong input. Please enter a number ...')
    
    ixmin = ixmin[iy]
    ixmax = ixmax[iy]
    L = ixmax - ixmin
    
    ys = yopts[ixmin:ixmax]
    
    # Define and plot the RCM-E fields for the chosen slice
    x_array = xin[ixmin:ixmax]
    # K_array = k[ixmin:ixmax]
    Vm_array = vm[ixmin:ixmax]
    # Cond_array = pedlam[ixmin:ixmax]
    
    Vgrad_array = np.zeros(L)
    for i in range(L):
        if(i > 0 and i < L - 1):
            Vgrad_array[i] = (Vm_array[i+1] - Vm_array[i-1])/(x_array[i+1] - x_array[i-1])
    
    # def K(x):
    #     y = np.interp(x,x_array,K_array)
    #     return y
    def Vgrad(x):
        y = np.interp(x,x_array,Vgrad_array)
        return y
    # def Cond(x):
    #     y = np.interp(x,x_array,Cond_array)
    #     return y
    
    # plt.figure(4); plt.clf()
    # plt.plot(x_array,K_array)
    # plt.grid(); plt.ylabel(r'$K(x)$'); plt.xlabel(r'$x$')
    # plt.title(r'Entropy at $T = $' + str(T) + ' and $Y = $' + str(ys))
    # plt.show()
    
    plt.figure(5); plt.clf()
    plt.plot(x_array,Vgrad_array)
    plt.grid(); plt.ylabel(r'$\nabla V^{-2/3} (x)$'); plt.xlabel(r'$x_i$')
    plt.title(r'Vgrad Term at $T = $' + str(T) + ' and $Y = $' + str(ys))
    plt.show()
    
    # plt.figure(6); plt.clf()
    # plt.plot(x_array,Cond_array)
    # plt.grid(); plt.ylabel(r'$\Sigma(x)$'); plt.xlabel(r'$x$')
    # plt.title(r'Pedersen Conductivity at $T = $' + str(T) + ' and $Y = $' + str(ys))
    # plt.show()
    
    # Plot alongside analytical profiles and fix parameters
    # Flux tube volume (fixes Scos = S * cos(zeta))
    
    Scos = -3/2*Vgrad(xspecial)
    Vgradanalytic = -2/3*Scos * np.ones(L)
    
    plt.figure(8); plt.clf()
    plt.plot(x_array,Vgrad_array,'b-',label = 'RCM-E')
    plt.plot(x_array,Vgradanalytic,'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\nabla V^{-2/3} (x)$'); plt.xlabel(r'$x_i$')
    plt.title('FTV Gradient Profile')
    plt.show()
    
    print('Note: These values only need to match at the X-slice of interest.')
    
    return Scos