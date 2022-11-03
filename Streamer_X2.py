import numpy as np
import matplotlib.pyplot as plt

def Streamer_X(T):
    
    # T = str(10)
    # Load in tab delimited file
    filein = 'Streamer_GRL2014_T='+T+'.dat'
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
    
    nbadL = 61
    nx = int(np.sqrt(xion.size)) - nbadL
    # higher nx causes problems because of chopped values in RCM
    
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
        plt.plot(yion[i,:],k[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('k'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(2); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],ftv[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(3); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],pedlam[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('pedlam'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(4); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],np.sqrt(brat[i,:]),label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('Sqrt brat'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(4); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],edistrat[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('Edistrat'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    # Choose X-slice
    try:
        ix = int(input('Enter your X-slice: ')) - 1
    except:
        print('Wrong input. Please enter a number ...')
    
    xs = xion[ix,0]
    
    # Define and plot the RCM-E fields for the chosen slice
    y_array = yion[ix,:]
    K_array = k[ix,:]
    Cond_array = pedlam[ix,:]
    Bmin_array = bmin[ix,:]
    
    def K(y):
        x = np.interp(y,y_array,K_array)
        return x
    def Cond(y):
        x = np.interp(y,y_array,Cond_array)
        return x
    
    plt.figure(6); plt.clf()
    plt.plot(y_array,K_array)
    plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Entropy at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    plt.figure(7); plt.clf()
    plt.plot(y_array,Cond_array)
    plt.grid(); plt.ylabel(r'$\Sigma(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Pedersen Conductivity at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    plt.figure(8); plt.clf()
    plt.plot(y_array,Bmin_array)
    plt.grid(); plt.ylabel(r'$Bmin(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Bmin at $T = $' + str(T) + ' and $X = $' + str(xs))
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
    L = 1000
    
    alpha = Kmin/K0
    eta = (1 - alpha)/(1 + alpha)
    
    r_array = np.linspace(r0,rf,L)
    
    def KR(r):
        x = (1 - eta*np.cos(np.pi*r))/(1 + eta)
        return x
    
    def KL(r):
        x = (1 - eta*np.cos(np.pi*r/r0))/(1 + eta)
        return x
    
    def Kanalytic(r):
        x = K0 * np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: KL(r), lambda r: KR(r), lambda r: 1])
        return x
    
    plt.figure(9); plt.clf()
    plt.plot(rnew,Knew,'b-',label = 'RCM-E')
    plt.plot(r_array,Kanalytic(r_array),'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$K(r)$'); plt.xlabel(r'$r$')
    plt.title(r'Entropy Profile at $x_s$ = ' + str(xs))
    plt.show()
    
    # Pedersen conductivity (fixes sigma)
    plt.figure(10); plt.clf()
    plt.plot(r_array,Cond_array,'b-')
    plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
    plt.title(r'Conductivity Profile at $x_s$ = ' + str(xs))
    plt.show()
    
    # Sigma0 = np.min(Cond_array)
    
    try:
        Sigma0 = float(input('Please enter the conductivity value corresponding to the flat region: '))
    except:
        print('Wrong input. Please enter a number ...')
    
    Sigmamax = np.max(Cond_array)
    sigma = (Sigmamax - Sigma0)/Sigma0
    
    try:
        beta = float(input('Please enter a guess for beta (larger values are thinner bumps): '))
    except:
        print('Wrong input. Please enter a number ...')
    
    def Cond(r):
        x = 1 + sigma * (np.sin(np.pi*r))**beta
        return x
    
    def Condanalytic(r):
        x = Sigma0 * np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: Cond(r), lambda r: 1])
        return x
    
    plt.figure(11); plt.clf()
    plt.plot(r_array,Cond_array,'b-',label = 'RCM-E')
    plt.plot(r_array,Condanalytic(r_array),'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
    plt.title(r'Conductivity Profile at $x_s$ = ' + str(xs))
    plt.show()

    return xs, d, K0, xi, eta, Sigma0, sigma, beta

# if you select wrong, re-prompt
# Automate 1 or both sides of y-boundary determination.
# Reprompt for beta
# Automate Sigma0
# Legends