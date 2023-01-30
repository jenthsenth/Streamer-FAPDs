#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import Solvers as s
import Streamer_X2 as sx
import Streamer_Y2 as sy
from textwrap import wrap

# Pick a specific streamer time and slices to determine parameters

try:
    T = input('Enter the moment of time you want to examine: ')
except:
    print('Wrong input. Please enter a number ...')

xs, d, K0, xi, eta, Sigma0, sigma, dyg, dm, Birk_array, RCMpot_array = sx.Streamer_X(T)
Scos, ys = sy.Streamer_Y(T,xs)

print('(xs,ys) = (' + str(xs) + ',' + str(ys) + ')')
print('The streamer width d = ' + str(d) + 'R_E')
print('The bubble width dmag = ' + str(dm) + 'R_E')
print('The entropy magnitude K0 = ' + str(K0) + ' nPa (R_E/nT)^5/3')
print('The asymmetry factor xi = ' + str(xi))
print('The factor eta = ' + str(eta))
print('The background Pedersen conductance is Sigma0 = ' + str(Sigma0) + ' S')
print('The conductivity enhancement factor is sigma = ' + str(sigma))
print('The FTV factor Scos(zeta) = ' + str(Scos) + ' 1/nT (R_E/nT)^(-5/3)')

# Define mapping ratios

drat = d/dm
ratio = drat

# Convert distance units from RE to m

RE = 6371009
d = d*RE
dyg = dyg*RE

# Define dimensionful scaling factors

Phi0 = (eta/(eta+1))*np.pi*K0*Scos*d/Sigma0
R0 = (d**2)/Sigma0
J0 = Phi0/R0
E0amp = Phi0/dyg

# FACs should be ~< 25 \muA/m^2
# PCP ~ 50 kV
# FAPD should be ~ 1-10 kV
# Magnetospheric field should be ~ 1-10 mV/m
# Ionospheric field should be ~ 10-50 mV/m

# Bubble width ~ 1-3 R_E
# Streamer width ~ 10-100 km (??)

# B_min ~< 10 nT
# v_BBF ~ 300-400 km/s
# v_str ~ 1 km/s


# automate slice selection to match to RCM potential drop

# solve phi equation first, make contour plot, determine direction and magnitude of field?

print('The dimensions for potential are',Phi0/1000,'kV.')
print('The dimensions for ionospheric electric field are',E0amp*1000,'mV / m.')
print('The dimensions for current density are',J0*(10**6),'\mu A / m^2.')
print('The dimensions for resistivity are',R0/10**6,'M \Omega m^2.')

# Define geometry
r_xi = (xi-1)/xi
r0 = r_xi
rf = 1
dr = (rf-r0)/1000

r_array = np.linspace(r0,rf,1000)

# Import the specific effective potential RHS to the RKA algorithm

def rhs(state,r):
    x = -state[0]*Veff(r) + Src(r)
    return x

# Create background profiles and potential

def f(r):
    x = np.pi*sigma * np.cos(np.pi*r) / (1 + sigma * np.sin(np.pi*r))
    return x

def Veff(r):
    x = np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < rf), r >= rf], [lambda r: 0, lambda r: f(r), lambda r: 0])
    return x

# plt.figure(1); plt.clf()
# plt.plot(r_array,Veff(r_array),'b-')
# plt.ylabel('V(r)'); plt.xlabel('r')
# plt.title('Potential')
# plt.grid('on')

def h(r):
    x = - xi * np.sin(np.pi*r) / (1 + sigma * np.sin(np.pi*r))
    return x

def g(r):
    x = - xi/r0 * np.sin(np.pi*r/r0)
    return x

def Src(r):
    x = np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < rf), r >= 1], [lambda r: 0, lambda r: g(r), lambda r: h(r), lambda r: 0])
    return x

# plt.figure(2); plt.clf()
# plt.plot(r_array,Src(r_array),'b-')
# plt.ylabel('Src(r)'); plt.xlabel('r')
# plt.title('Source Term')
# plt.grid('on')

def H(r):
    x = - np.sin(np.pi*r)/xi
    return x

def J(r):
    x = np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 0, lambda r: g(r)/xi**2, lambda r: H(r), lambda r: 0])
    return x

# Plot FACs

plt.figure(4); plt.clf()
plt.plot(r_array,J0*J(r_array)*10**6,'r-',label = 'Analytic')
plt.plot(r_array,2*Birk_array,'b-',label = 'RCM-E')
plt.legend(); plt.grid('on'); plt.ylabel(r'$J_\parallel(r) \ \ (\ \mu A \, / m^2)$'); plt.xlabel('r')
plt.title('Field-Aligned Current')



# Use adaptive Runge-Kutta to solve for ionospheric electric field. Prime ' is r-derivative:

# Solver for ionospheric electric field
E0 = 0.0
r = r0
rplot = np.array([])
drplot = np.array([])
Eplot = np.array([])
state = np.array([E0])
adaptErr = 1.0e-12

while r < rf:
    [state,r,dr] = s.rka(state,r,dr,adaptErr,rhs)
    # store plots
    rplot = np.append(rplot,r)
    drplot = np.append(drplot,dr)
    Eplot = np.append(Eplot,state[0])

plt.figure(5); plt.clf()
plt.plot(rplot,E0amp*Eplot*1000,'b-')
plt.ylabel(r'$E_i(r) \ (mV/m)$'); plt.xlabel('r')
plt.title('Ionospheric Electric Field at T = ' + str(T))
plt.grid('on')

# Compute FAPD

R = (1/R0)*(10**6)*np.float(input('Enter the resistivity (in M\Omega m^2): '))

Jplot = J(rplot)

def Grad_J(r):
    x = np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 0, lambda r: j(r), lambda r: k(r), lambda r: 0])
    return x

def j(r):
    x = - (np.pi/r0**2) * np.cos(np.pi*r/r0) / xi
    return x

def k(r):
    x = - np.pi * np.cos(np.pi*r) / xi
    return x

GradJplot = Grad_J(rplot)*np.where(Jplot < 0, 1, 0)

plt.figure(6); plt.clf()
plt.plot(rplot,E0amp*R*GradJplot/1000,'b-')
plt.ylabel(r'$\Delta \Phi(r) \, (kV)$'); plt.xlabel('r')
plt.title("\n".join(wrap('Field-Aligned Potential Drop at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
plt.grid('on')

# Compute magnetospheric electric field

Emplot = Eplot - R*GradJplot

plt.figure(7); plt.clf()
plt.plot(rplot,ratio*E0amp*Emplot*1000,'r-')
plt.ylabel(r'$E_m(r) \ (mV/m)$'); plt.xlabel('r')
plt.title("\n".join(wrap('Magnetospheric Electric Field at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
plt.grid('on')

# Compute potentials

L = rplot.size
potential = np.zeros(L)
potentialm = np.zeros(L)
potentialm1 = np.zeros(L)
for i in range(L-1):
    dr = rplot[i+1]-rplot[i]
    potential[i+1] = potential[i]-dr*(Eplot[i]+Eplot[i+1])/2
    potentialm1[i+1] = potentialm1[i]-dr*(Emplot[i]+Emplot[i+1])/2

potentialm = potential + R*Jplot*np.where(Jplot < 0, 1, 0)

plt.figure(9); plt.clf()
plt.plot(rplot,Phi0*potential/1000,'b-')
plt.plot(rplot,Phi0*potentialm/1000,'r-')
plt.plot(r_array,RCMpot_array/1000 - RCMpot_array[0]/1000,'g-')
plt.ylabel(r'$\Phi(r) \, (kV)$'); plt.xlabel('r')
plt.legend(['potential', 'potentialm','potentialRCM'])
plt.title("\n".join(wrap('Magnetospheric and Ionospheric Potentials at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
plt.grid('on')



# SANITY CHECKS

# Jplot = J(rplot)
# Jplot[Jplot > 0] = 0

# # plt.figure(8); plt.clf()
# # plt.plot(rplot,Phi0*R*Jplot/1000,'b-')
# # plt.ylabel(r'$\Delta \Phi(r) \, (kV)$'); plt.xlabel('r')
# # plt.title("\n".join(wrap('Field-Aligned Potential Drop at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
# # plt.grid('on')

# DJplot = np.zeros(L)
# for i in range(L - 1):
#     if(i > 0 and i < L - 2):
#         DJplot[i] = (-Jplot[i+2] + 8*Jplot[i+1] - 8*Jplot[i-1] + Jplot[i-2])/(6*(rplot[i+1] - rplot[i-1]))

# plt.figure(9); plt.clf()
# plt.plot(rplot,E0amp*R*DJplot/1000,'b-')
# plt.title('DJplot')
# plt.grid('on')



# Emplot = Eplot - R*DJplot


# plt.figure(7); plt.clf()
# plt.plot(rplot,ratio*E0amp*Emplot1*1000,'b-')
# plt.plot(rplot,ratio*E0amp*Emplot*1000,'r-')
# plt.ylabel(r'$E_m(r) \ (mV/m)$'); plt.xlabel('r')
# plt.title("\n".join(wrap('Magnetospheric Electric Field at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
# plt.grid('on')

# DEplot = np.zeros(L)
# for i in range(L - 1):
#     if(i > 0 and i < L - 2):
#         DEplot[i] = (-Eplot[i+2] + 8*Eplot[i+1] - 8*Eplot[i-1] + Eplot[i-2])/(6*(rplot[i+1] - rplot[i-1]))

# def DEplot(r):
#     x = np.interp(r,rplot,DE_array)
#     return x

# Zeroplot = DEplot + Veff(rplot)*Eplot - Src(rplot)

# plt.figure(11); plt.clf()
# plt.plot(rplot,Zeroplot,'b-')
# plt.title('Equation Satisfied if Near Zero')
# plt.grid('on')

# plt.figure(10); plt.clf()
# plt.plot(rplot,E0amp*R*GradJplot/1000,'b-')
# # plt.plot(rplot,E0amp*R*DJplot/1000,'r-')
# plt.ylabel(r'$\Delta \Phi(r) \, (kV)$'); plt.xlabel('r')
# plt.title("\n".join(wrap('Field-Aligned Potential Drop at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
# plt.grid('on')


# plt.figure(9); plt.clf()
# plt.plot(rplot,Phi0*potential/1000,'b-')
# # plt.plot(rplot,Phi0*potentialm1/1000,'g-')
# plt.plot(rplot,Phi0*potentialm/1000,'r-')
# plt.ylabel(r'$\Phi(r) \, (kV)$'); plt.xlabel('r')
# plt.legend(['potential', 'potentialm'])
# plt.title("\n".join(wrap('Magnetospheric and Ionospheric Potentials at T = ' + str(T) + r' for Resistivity $R$ = ' + str(R*R0/10**6) + r' $M\Omega \, m^2$',60)))
# plt.grid('on')

# Use sys to make folder architecture better