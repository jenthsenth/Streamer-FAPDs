import numpy as np
import matplotlib.pyplot as plt

def Streamer_X(T):
    
    # Load in tab delimited file
    filein = 'T'+T+'xi.txt'
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
    
    # Find unique values of x
    xopts = xin[0]
    xlast = xopts
    
    # Count unique X-slices and determine corresponding y-indices
    iy = np.array(0)
    
    for i in range(xin.size):
        if(xin[i]!=xlast):
            iy = np.append(iy,i)
            xopts = np.append(xopts,xin[i])
            xlast = xin[i]
    
    iy = np.append(iy,xin.size)
    ix = xopts.size
    
    print('There are ',ix,' unique x values in this dataset')
    print('x = ',xopts)
    
    # Plot the RCM-E fields for all slices at the given time step
    iymin = iy[0:ix]
    iymax = iy[1:]
    
    plt.figure(1); plt.clf()
    for i in range(xopts.size):
        plt.plot(yin[iymin[i]:iymax[i]],k[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
    plt.legend(); plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(2); plt.clf()
    for i in range(xopts.size):
        plt.plot(yin[iymin[i]:iymax[i]],ftv[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
    plt.legend(); plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(3); plt.clf()
    for i in range(xopts.size):
        plt.plot(yin[iymin[i]:iymax[i]],pedlam[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
    plt.legend(); plt.grid(); plt.ylabel('pedlam'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    # Choose X-slice
    try:
        ix = int(input('Enter your X-slice: ')) - 1
    except:
        print('Wrong input. Please enter a number ...')
    
    xspecial = xopts[ix]
    
    iymin = iymin[ix]
    iymax = iymax[ix]
    xs = xopts[iymin:iymax]
    
    # Define and plot the RCM-E fields for the chosen slice
    y_array = yin[iymin:iymax]
    K_array = k[iymin:iymax]
    Cond_array = pedlam[iymin:iymax]
    
    def K(y):
        x = np.interp(y,y_array,K_array)
        return x
    def Cond(y):
        x = np.interp(y,y_array,Cond_array)
        return x
    
    plt.figure(4); plt.clf()
    plt.plot(y_array,K_array)
    plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Entropy at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    plt.figure(6); plt.clf()
    plt.plot(y_array,Cond_array)
    plt.grid(); plt.ylabel(r'$\Sigma(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Pedersen Conductivity at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    # Shift y-boundaries
    y0new = float(input('Enter a new y_0: '))
    yfnew = float(input('Enter a new y_f: '))
    
    ynew = np.linspace(y0new,yfnew,1000)
    d = yfnew - y0new
    
    Knew = K(ynew)
    Cond_array = Cond(ynew)
    
    # Rescale and shift y-axis to create r-axis
    Kmin = np.min(Knew)
    imin = np.argmin(Knew)
    ymin = ynew[imin]
    rescale = ynew[-1] - ymin
    rnew = ynew/rescale
    rmin = ymin/rescale
    rmin_array = rmin*np.ones(1000)
    rnew = rnew - rmin_array
    
    # Plot alongside analytical profiles and fix parameters
    # Entropy (fixes K0, xi, eta, d)
    K0 = np.max(Knew)
    
    r0 = rnew[0]
    xi = 1/(1 - r0)
    rf = 1
    dr = (rf-r0)/1000
    
    alpha = Kmin/K0
    eta = (1 - alpha)/(1 + alpha)
    
    r_array = np.arange(r0,rf,dr)
    
    def KR(r):
        x = (1 - eta*np.cos(np.pi*r))/(1 + eta)
        return x
    
    def KL(r):
        x = (1 - eta*np.cos(np.pi*r/r0))/(1 + eta)
        return x
    
    def Kanalytic(r):
        x = K0 * np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: KL(r), lambda r: KR(r), lambda r: 1])
        return x
    
    plt.figure(7); plt.clf()
    plt.plot(rnew,Knew,'b-',label = 'RCM-E')
    plt.plot(r_array,Kanalytic(r_array),'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$K(r)$'); plt.xlabel(r'$r$')
    plt.title('Entropy Profile')
    plt.show()
    
    # Pedersen conductivity (fixes sigma)
    plt.figure(10); plt.clf()
    plt.plot(r_array,Cond_array,'b-')
    plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
    plt.title('Conductivity Profile')
    plt.show()
    
    try:
        Sigma0 = float(input('Please enter the conductivity value corresponding to the flat region: '))
    except:
        print('Wrong input. Please enter a number ...')
    
    Sigmamax = np.max(Cond_array)
    sigma = (Sigmamax - Sigma0)/Sigma0
    
    beta = float(input('Please enter a guess for beta (larger values are thinner bumps): '))
    
    def Cond(r):
        x = 1 + sigma * (np.sin(np.pi*r))**beta
        return x
    
    def Condanalytic(r):
        x = Sigma0 * np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: Cond(r), lambda r: 1])
        return x
    
    plt.figure(10); plt.clf()
    plt.plot(r_array,Cond_array,'b-',label = 'RCM-E')
    plt.plot(r_array,Condanalytic(r_array),'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
    plt.title('Conductivity Profile')
    plt.show()
    
    return xspecial, d, K0, xi, eta, Sigma0, sigma, beta

# if you select wrong, re-prompt
# Automate 1 or both sides of y-boundary determination.
# Reprompt for beta
# Automate Sigma0
# Legends