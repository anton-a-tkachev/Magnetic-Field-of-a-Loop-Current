#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 19:02:25 2020

@author: atkachev
"""

import numpy as np
from magnetic import MF
import matplotlib.pyplot as plt
from scipy.integrate import odeint, RK45, solve_ivp

eta = 175882001077.216
e = 1.602176634E-19
me = 9.1093837015E-31
E0 = 510998.949996164
U = 30E3
c = 299792458
gamma = 1 + U/E0
beta = np.sqrt(1 - 1/gamma**2)

COILS = [
    {
        'Label' :   'EMITTER COIL',
        'Z0'    :   0.200,
        'ID'    :   0.550,
        'OD'    :   0.600,
        'DZ'    :   0.050,
        'W'     :   200,
        'I'     :   150,
        'NR'    :   2,
        'NZ'    :   5
    },
    {
        'Label' :   'BEAMLINE SOLENOID',
        'Z0'    :   0.800,
        'ID'    :   0.250,
        'OD'    :   0.300,
        'DZ'    :   0.250,
        'W'     :   400,
        'I'     :   25, 
        'NR'    :   2,
        'NZ'    :   10
    },
    {
        'Label' :   'NEW MIRROR PLUG',
        'Z0'    :   2.500,
        'ID'    :   0.800,
        'OD'    :   0.900,
        'DZ'    :   0.200,
        'W'     :   1000,
        'I'     :   1400, 
        'NR'    :   3,
        'NZ'    :   10
    }
]

def ABC(R, Z, COILS):
    PP = np.array(((1,1,0),))
    for COIL in COILS:
        II = COIL['I']*COIL['W']/COIL['NR']/COIL['NZ']*np.ones((COIL['NR']*COIL['NZ'], 1))
        ZZ = np.linspace(-0.5*COIL['DZ'], 0.5*COIL['DZ'], COIL['NZ']) + COIL['Z0']
        AA = np.linspace(0.5*COIL['ID'], 0.5*COIL['OD'], COIL['NR'])
        PP = np.concatenate((PP, np.concatenate((np.transpose([np.tile(AA, len(ZZ)), np.repeat(ZZ, len(AA))]), II), axis=1)), axis=0)
    return np.sum([MF(R, Z - pp[1], pp[0], pp[2]) for pp in PP[1:,:]], axis=0)

zlist = np.linspace(0.0, 3.0, 300)
rlist = np.linspace(0.001, 0.1, 100)
R, Z = np.meshgrid(rlist, zlist)
F = ABC(R, Z, COILS)[0]
fig,ax=plt.subplots(1,1)
lvls = np.linspace(0.0, 0.002, 41)
cp = ax.contour(Z*100, R*100, F, levels=lvls, linewidths=0.5)
#fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Azimuthal vector-potential (Wb/m)')
ax.set_xlabel('z (cm)')
ax.set_ylabel('r (cm)')
plt.show()

p0 = 0.0
R0 = 0.010
Z0 = 0.000
VR0 = 0.000
VZ0 = beta*c
P0 = p0/e + ABC(R0,Z0,COILS)[0]
x0 = [R0, Z0, VR0, VZ0]

def fun(t, x):
    [A, Br, Bz] = ABC(x[0], x[1], COILS)
    return [x[2], x[3], (eta/gamma)**2*(P0 - A)*Bz, (eta/gamma)**2*(P0 - A)*Br]     # xdot

solution = solve_ivp(fun, (0, 1E-8), x0, method='RK45', max_step=1e-12)
print(solution['message'])
sol = solution['y']
R = 100*sol[0,:]
Z = 100*sol[1,:]
plt.plot(Z, R, color='red', linewidth=0.5)
plt.xlim(0,300)
plt.ylim(0,10)
plt.show()

# dt = 7.111E-12
# p0 = 0.0
# R0 = 0.010
# Z0 = 0.000
# X0 = np.array([[R0, Z0],])
# X1 = X0 + np.array([[0, beta*c*dt],])
# P0 = p0/e + ABC(R0,Z0,COILS)[0]
# XX = np.concatenate((X0, X1), axis=0)

# for i in range(2500):
#     R, Z = X1[0,0], X1[0,1]
#     [A, Br, Bz] = ABC(R, Z, COILS)
#     a = (eta/gamma)**2*(P0 - A)*np.array([[Bz, Br],])
#     X2 = 2*X1 - X0 + a*dt**2
#     XX = np.concatenate((XX,X2), axis=0)
#     X0, X1 = X1, X2