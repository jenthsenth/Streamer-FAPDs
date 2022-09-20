import numpy as np
import matplotlib.pyplot as plt

# def Streamer_X(T):

# Load in tab delimited file
filein = 'Streamer File.dat'
data = np.loadtxt(filein,skiprows=65)

# Read in field arrays
xmag = data[:,6]
ymag = data[:,7]
vm = data[:,8]
ftv = data[:,9]
v = data[:,10]
birk = data[:,12]
pedlam = data[:,12]
hall = data[:,16]
eflux = data[:,17]
N = data[:,19]
T = data[:,20]
P = data[:,21]
k = data[:,22]
Ne = data[:,23]
Te = data[:,24]
Pe = data[:,25]
ke = data[:,26]
Ni = data[:,27]
Ti = data[:,28]
Pi = data[:,29]
ki = data[:,30]
xion = data[:,43]
yion = data[:,44]

n = 10

xmag = xmag[0::n]
ymag = ymag[0::n]
vm = vm[0::n]
ftv = ftv[0::n]
v = v[0::n]
birk = birk[0::n]
pedlam = pedlam[0::n]
hall = hall[0::n]
eflux = eflux[0::n]
N = N[0::n]
T = T[0::n]
P = P[0::n]
k = k[0::n]
Ne = Ne[0::n]
Te = Te[0::n]
Pe = Pe[0::n]
ke = ke[0::n]
Ni = Ni[0::n]
Ti = Ti[0::n]
Pi = Pi[0::n]
ki = ki[0::n]
xion = xion[0::n]
yion = yion[0::n]

# Find unique values of x
xopts = xion[0]
xlast = xopts

# Count unique X-slices and determine corresponding y-indices
iy = np.array(0)

for i in range(xion.size):
    if(xion[i]!=xlast):
        iy = np.append(iy,i)
        xopts = np.append(xopts,xion[i])
        xlast = xion[i]

iy = np.append(iy,xion.size)
ix = xopts.size

print('There are ',ix,' unique x values in this dataset')
print('x = ',xopts)

# Plot the RCM-E fields for all slices at the given time step
iymin = iy[0:ix]
iymax = iy[1:]

plt.figure(1); plt.clf()
for i in range(xopts.size):
    plt.plot(yion[iymin[i]:iymax[i]],k[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
plt.legend(); plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y_i$')
# plt.title(filein)
plt.show()

plt.figure(2); plt.clf()
for i in range(xopts.size):
    plt.plot(yion[iymin[i]:iymax[i]],ftv[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
plt.legend(); plt.grid(); plt.ylabel('ftv'); plt.xlabel(r'$y_i$')
plt.title(filein)
plt.show()

plt.figure(3); plt.clf()
for i in range(xopts.size):
    plt.plot(yion[iymin[i]:iymax[i]],pedlam[iymin[i]:iymax[i]],label='x = ' + str(xopts[i]))
plt.legend(); plt.grid(); plt.ylabel('pedlam'); plt.xlabel(r'$y_i$')
plt.title(filein)
plt.show()

# # Choose X-slice
# try:
#     ix = int(input('Enter your X-slice: ')) - 1
# except:
#     print('Wrong input. Please enter a number ...')

# xspecial = xopts[ix]

# iymin = iymin[ix]
# iymax = iymax[ix]
# xs = xopts[iymin:iymax]

# # Define and plot the RCM-E fields for the chosen slice
# y_array = yion[iymin:iymax]
# K_array = k[iymin:iymax]
# Cond_array = pedlam[iymin:iymax]

# def K(y):
#     x = np.interp(y,y_array,K_array)
#     return x
# def Cond(y):
#     x = np.interp(y,y_array,Cond_array)
#     return x

# plt.figure(4); plt.clf()
# plt.plot(y_array,K_array)
# plt.grid(); plt.ylabel(r'$K(y)$'); plt.xlabel(r'$y$')
# plt.title(r'Entropy at $T = $' + str(T) + ' and $X = $' + str(xs))
# plt.show()

# plt.figure(6); plt.clf()
# plt.plot(y_array,Cond_array)
# plt.grid(); plt.ylabel(r'$\Sigma(y)$'); plt.xlabel(r'$y$')
# plt.title(r'Pedersen Conductivity at $T = $' + str(T) + ' and $X = $' + str(xs))
# plt.show()

# # Shift y-boundaries
# y0new = float(input('Enter a new y_0: '))
# yfnew = float(input('Enter a new y_f: '))

# ynew = np.linspace(y0new,yfnew,1000)
# d = yfnew - y0new

# Knew = K(ynew)
# Cond_array = Cond(ynew)

# # Rescale and shift y-axis to create r-axis
# Kmin = np.min(Knew)
# imin = np.argmin(Knew)
# ymin = ynew[imin]
# rescale = ynew[-1] - ymin
# rnew = ynew/rescale
# rmin = ymin/rescale
# rmin_array = rmin*np.ones(1000)
# rnew = rnew - rmin_array

# # Plot alongside analytical profiles and fix parameters
# # Entropy (fixes K0, xi, eta, d)
# K0 = np.max(Knew)

# r0 = rnew[0]
# xi = 1/(1 - r0)
# rf = 1
# dr = (rf-r0)/1000

# alpha = Kmin/K0
# eta = (1 - alpha)/(1 + alpha)

# r_array = np.arange(r0,rf,dr)

# def KR(r):
#     x = (1 - eta*np.cos(np.pi*r))/(1 + eta)
#     return x

# def KL(r):
#     x = (1 - eta*np.cos(np.pi*r/r0))/(1 + eta)
#     return x

# def Kanalytic(r):
#     x = K0 * np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: KL(r), lambda r: KR(r), lambda r: 1])
#     return x

# plt.figure(7); plt.clf()
# plt.plot(rnew,Knew,'b-',label = 'RCM-E')
# plt.plot(r_array,Kanalytic(r_array),'r-', label = 'Analytic')
# plt.legend(); plt.grid('on'); plt.ylabel(r'$K(r)$'); plt.xlabel(r'$r$')
# plt.title('Entropy Profile')
# plt.show()

# # Pedersen conductivity (fixes sigma)
# plt.figure(10); plt.clf()
# plt.plot(r_array,Cond_array,'b-')
# plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
# plt.title('Conductivity Profile')
# plt.show()

# try:
#     Sigma0 = float(input('Please enter the conductivity value corresponding to the flat region: '))
# except:
#     print('Wrong input. Please enter a number ...')

# Sigmamax = np.max(Cond_array)
# sigma = (Sigmamax - Sigma0)/Sigma0

# beta = float(input('Please enter a guess for beta (larger values are thinner bumps): '))

# def Cond(r):
#     x = 1 + sigma * (np.sin(np.pi*r))**beta
#     return x

# def Condanalytic(r):
#     x = Sigma0 * np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 1, lambda r: Cond(r), lambda r: 1])
#     return x

# plt.figure(10); plt.clf()
# plt.plot(r_array,Cond_array,'b-',label = 'RCM-E')
# plt.plot(r_array,Condanalytic(r_array),'r-', label = 'Analytic')
# plt.legend(); plt.grid('on'); plt.ylabel(r'$\Sigma(r)$'); plt.xlabel(r'$r$')
# plt.title('Conductivity Profile')
# plt.show()

# # return xspecial, d, K0, xi, eta, Sigma0, sigma, beta

# # if you select wrong, re-prompt
# # Automate 1 or both sides of y-boundary determination.
# # Reprompt for beta
# # Automate Sigma0
# # Legends