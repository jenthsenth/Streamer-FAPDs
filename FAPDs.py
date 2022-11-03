#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import Solvers as s
import Streamer_X2 as sx
import Streamer_Y2 as sy
from textwrap import wrap

# Pick a specific streamer time and scenario
try:
    T = input('Enter the moment of time you want to examine: ')
except:
    print('Wrong input. Please enter a number ...')

xs, d, K0, xi, eta, Sigma0, sigma, beta = sx.Streamer_X(T)
Scos, ys = sy.Streamer_Y(T,xs)

print('(xs,ys) = (' + str(xs) + ',' + str(ys) + ')')
print('The streamer width d = ' + str(d) + 'R_E')
print('The entropy magnitude K0 = ' + str(K0) + ' nPa (R_E/nT)^5/3')
print('The asymmetry factor xi = ' + str(xi))
print('The factor eta = ' + str(eta))
print('The background Pedersen conductance is Sigma0 = ' + str(Sigma0) + ' S')
print('The conductivity enhancement factor is sigma = ' + str(sigma))
print('The conductivity width factor is beta = ' + str(beta))
print('The FTV factor Scos(zeta) = ' + str(Scos) + ' R_E/nT')

RE = 6371009
Phi0 = (eta/(eta+1))*np.pi*K0*Scos*d*RE/Sigma0
E0amp = Phi0/RE*1000
J0 = Phi0*Sigma0/(d**2)/(RE**2)*(10**6)
R0 = (d**2)*(RE**2)/Sigma0/10**6

ratio = 0.1

print('The dimensions for ionospheric electric field are',E0amp,'mV / m.')
print('The dimensions for current density are',J0,'\mu A / m^2.')
print('The dimensions for resistivity are',R0,'M \Omega m^2.')

# Import the specific effective potential RHS to the RKA algorithm
def rhs(state,r):
    x = -state[0]*Veff(r) + Src(r)
    return x

# Define radial domain & other parameters
adaptErr = 1.0e-5
r_xi = (xi-1)/xi
r0 = r_xi
rf = 1
dr = (rf-r0)/1000

r_array = np.arange(r0,rf,dr)

# Create background profiles and potential
def f(r):
    x = np.pi*beta*sigma * ((np.sin(np.pi*r))**(beta - 1)) * np.cos(np.pi*r) / (1 + sigma * (np.sin(np.pi*r))**beta)
    return x

def Veff(r):
    x = np.piecewise(r, [r < 0, np.logical_and(r >= 0, r < rf), r >= rf], [lambda r: 0, lambda r: f(r), lambda r: 0])
    return x

# plt.figure(2); plt.clf()
# plt.plot(r_array,Veff(r_array),'b-')
# plt.ylabel('V(r)'); plt.xlabel('r')
# plt.title('Potential')
# plt.grid('on')

def h(r):
    x = - xi * np.sin(np.pi*r) / (1 + sigma * (np.sin(np.pi*r))**beta)
    return x

def g(r):
    x = - xi/r0 * np.sin(np.pi*r/r0)
    return x

def Src(r):
    x = np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 0, lambda r: g(r), lambda r: h(r), lambda r: 0])
    return x

# plt.figure(3); plt.clf()
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

# plt.figure(4); plt.clf()
# plt.plot(r_array,J0*J(r_array),'b-')
# plt.ylabel(r'$J_\parallel(r) \ \ (\ \mu A \,/ m^2)$'); plt.xlabel('r')
# plt.title('Field-Aligned Current')
# plt.grid('on')

R = (1/R0)*np.float(input('Enter the resistivity (in M\Omega m^2): '))

# Perform RK4 algorithm to construct perturbation. Prime ' is r-derivative:
# u'   =  w
# w' = -Veff(r)*u

# Solver for ionospheric electric field
E0 = 0.0
r = r0
rplot = np.array([])
drplot = np.array([])
Eplot = np.array([])
state = np.array([E0])

while r < rf:
    [state,r,dr] = s.rka(state,r,dr,adaptErr,rhs)
    # store plots
    rplot = np.append(rplot,r)
    drplot = np.append(drplot,dr)
    Eplot = np.append(Eplot,state[0])

plt.figure(5); plt.clf()
plt.plot(rplot,E0amp*Eplot,'b-')
plt.ylabel(r'$E_i(r) \ (mV/m)$'); plt.xlabel('r')
plt.title('Ionospheric Electric Field at T = '+str(T))
plt.grid('on')

def Grad_J(r):
    x = np.piecewise(r, [r < r0, np.logical_and(r >= r0, r < 0), np.logical_and(r >= 0, r < 1), r >= 1], [lambda r: 0, lambda r: j(r), lambda r: k(r), lambda r: 0])
    return x

def j(r):
    x = - (np.pi/r0**2) * np.cos(np.pi*r/r0) / xi
    return x

def k(r):
    x = - np.pi * np.cos(np.pi*r) / xi
    return x

GradJplot = Grad_J(rplot)
GradJplot[GradJplot > 0] = 0

Emplot = Eplot + R*GradJplot

plt.figure(7); plt.clf()
plt.plot(rplot,ratio*E0amp*Emplot,'b-')
plt.ylabel(r'$E_m(r) \ (mV/m)$'); plt.xlabel('r')
plt.title("\n".join(wrap('Magnetospheric Electric Field at T = '+str(T)+r' for Resistivity $R$ = '+str(R*R0)+r' $M\Omega \, m^2$',60)))
plt.grid('on')

# plt.figure(8); plt.clf()
# plt.plot(rplot,R0*J0*R*Grad_J(rplot)/d,'b-')
# plt.ylabel(r'$\Delta \Phi(r)$'); plt.xlabel('r')
# plt.title('Field-Aligned Potential Drop')
# plt.grid('on')

# Use sys to make folder architecture better