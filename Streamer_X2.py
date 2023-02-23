import numpy as np
import matplotlib.pyplot as plt

def Streamer_X(T,Re):
    
    # Load in tab delimited file
    filein = 'Streamer_GRL2014_T='+T+'.dat'
    data = np.loadtxt(filein,skiprows=77)
    
    # filein = 'Streamer_WQ_T='+T+'.dat'
    # data = np.loadtxt(filein,skiprows=46)
    
    # Read in field arrays
    # xline = 44
    xline = 38
    
    xion = data[:,xline]
    nslices = 100
    ndim = xion.size
    
    xmag = data[:,6].reshape((int(ndim/nslices),nslices))
    ymag = data[:,7].reshape((int(ndim/nslices),nslices))
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
    yion = data[:,45].reshape((int(ndim/nslices),nslices))
    
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
    
    nbad = 61
    # nbad = 48
    nx = int(np.sqrt(xion.size)) - nbad
    # higher nx causes problems because of chopped values in RCM
    
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
    
    Robinson = (40*eavg_para/(16 + eavg_para**2))*np.sqrt(eflux_para)
    
    # Plot the RCM-E fields for all slices at the given time step
    plt.figure(1); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],k[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel('K'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(2); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],birk[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel(r'$J_\parallel \left[ \mu A/m^2 \right]$'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(3); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],pedpsi[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel(r'$\Sigma_{\phi \phi,diff} \left[ S \right]$'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(4); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],eflux[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel(r'$\bar{\Phi}_\parallel \left[ erg/cm^2 s \right]$'); plt.xlabel(r'$y_i$')
    plt.title(filein)
    plt.show()
    
    plt.figure(5); plt.clf()
    for i in range(nx):
        plt.plot(yion[i,:],Robinson[i,:],label='x = ' + str(xion[i,0]))
    plt.grid(); plt.ylabel(r'$\Sigma_{\phi \phi,mono} \left[ S \right]$'); plt.xlabel(r'$y_i$')
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
    ym_array = ymag[ix,:]
    K_array = k[ix,:]
    Cond0_array = pedpsi[ix,:]
    v_array = v[ix,:]
    Birk_array = birk[ix,:]
    # Eflux_array = eflux[ix,:]
    # Eavg_array = eavg[ix,:]
    # P_array = P[ix,:]
    # Ne_array = Ne[ix,:]
    # Tempe_array = Tempe[ix,:]
    # Efluxpara_array = eflux_para[ix,:]
    # Eavgpara_array = eavg_para[ix,:]
    Robinson_array = Robinson[ix,:]
    
    def ymag(y):
        x = np.interp(y,y_array,ym_array)
        return x
    def K(y):
        x = np.interp(y,y_array,K_array)
        return x
    def Cond0(y):
        x = np.interp(y,y_array,Cond0_array)
        return x
    def RCMpot(y):
        x = np.interp(y,y_array,v_array)
        return x
    def Birk(y):
        x = np.interp(y,y_array,Birk_array)
        return x
    # def Eflux(y):
    #     x = np.interp(y,y_array,Eflux_array)
    #     return x
    # def Eavg(y):
    #     x = np.interp(y,y_array,Eavg_array)
    #     return x
    # def P(y):
    #     x = np.interp(y,y_array,P_array)
    #     return x
    # def Ne(y):
    #     x = np.interp(y,y_array,Ne_array)
    #     return x
    # def Tempe(y):
    #     x = np.interp(y,y_array,Tempe_array)
    #     return x
    # def Eflux_para(y):
    #     x = np.interp(y,y_array,Efluxpara_array)
    #     return x
    # def Eavg_para(y):
    #     x = np.interp(y,y_array,Eavgpara_array)
    #     return x
    def CondEnhance(y):
        x = np.interp(y,y_array,Robinson_array)
        return x
    
    plt.figure(6); plt.clf()
    plt.plot(y_array,K_array)
    plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y$')
    plt.title(r'Entropy at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    plt.figure(7); plt.clf()
    plt.plot(y_array,Birk_array)
    plt.grid(); plt.ylabel(r'$J_\parallel(y) \, \left( \mu A / m^2 \right)$'); plt.xlabel(r'$y$')
    plt.title(r'Birkeland Currents at $T = $' + str(T) + ' and $X = $' + str(xs))
    plt.show()
    
    # Shift y-boundaries
    y0new = float(input('Enter a new y_0 (for LHS of entropy depletion): '))
    yfnew = float(input('Enter a new y_f (for RHS of entropy depletion): '))
    
    ynew = np.linspace(y0new,yfnew,1000)
    d = yfnew - y0new
    dm = ymag(yfnew) - ymag(y0new)
    
    Knew = K(ynew)
    Cond0_array = Cond0(ynew)
    Birk_array = Birk(ynew)
    RCMpot_array = RCMpot(ynew)
    CondEnhance_array = CondEnhance(ynew)
    Condtot_array = Cond0_array + CondEnhance_array
    
    # Rescale and shift y-axis to create r-axis
    Kmin = np.min(Knew)
    imin = np.argmin(Knew)
    ymin = ynew[imin]
    dyg = ynew[-1] - ymin
    rnew = ynew/dyg
    rmin = ymin/dyg
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
    
    plt.figure(8); plt.clf()
    plt.plot(rnew,Knew,'b-',label = 'RCM-E')
    plt.plot(r_array,Kanalytic(r_array),'r-',label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$K(r)$'); plt.xlabel(r'$r$')
    plt.title(r'Entropy Profile at $x_s$ = ' + str(xs))
    plt.show()
    
    # Pedersen conductivity (fixes sigma)
    
    # Sigma0 = np.average(Cond0_array)
    Sigma0 = 12
    
    # try:
    #     Sigma0 = float(input('Please enter the conductivity value corresponding to the flat region (or minimum conductivity): '))
    # except:
    #     print('Wrong input. Please enter a number ...')
    
    Sigmamax = np.max(Condtot_array)
    sigma = (Sigmamax - Sigma0)/Sigma0
    
    def Cond(r):
        x = 1 + sigma * (np.sin(np.pi*r))
        return x
    
    def Condanalytic(r):
        x = Sigma0 * np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: Cond(r), lambda r: 1])
        return x
    
    plt.figure(9); plt.clf()
    plt.plot(r_array,Condtot_array,'b-',label = 'RCM-E')
    plt.plot(r_array,Condanalytic(r_array),'r-', label = 'Analytic')
    plt.legend(); plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
    plt.title(r'Conductivity Profile at $x_s$ = ' + str(xs))
    plt.show()

    return xs, d, K0, xi, eta, Sigma0, sigma, dyg, dm, Birk_array, RCMpot_array

# if you select wrong, re-prompt
# Automate 1 or both sides of y-boundary determination.
# Automate Sigma0
# Legends