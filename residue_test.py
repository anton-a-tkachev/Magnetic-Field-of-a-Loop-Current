# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 17:30:04 2020

@author: atkachev
"""

from numpy import cos, linspace, array, sqrt
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
from scipy.special import ellipk, ellipe
from scipy.constants import pi, mu_0

dBz = lambda p,r,z,a: (a**2 - r*a*cos(p))/(r**2 + z**2 + a**2 - 2*r*a*cos(p))**(1.5)
Bz1 = lambda r,z,a,I: (mu_0*I/4/pi)*quad(dBz, 0, 2*pi, args=(r,z,a))[0]

dBr = lambda p,r,z,a: a*z*cos(p)/(r**2 + z**2 + a**2 - 2*r*a*cos(p))**(1.5)
Br1 = lambda r,z,a,I: (mu_0*I/4/pi)*quad(dBr, 0, 2*pi, args=(r,z,a))[0]

k = lambda r, z, a: 4*a*r/((a + r)**2 + z**2)
Bz2 = lambda r, z, a, i: mu_0*i/2/pi/sqrt((a + r)**2 + z**2)*(ellipk(k(r,z,a)) + (a**2 - r**2 - z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))
Br2 = lambda r, z, a, i: mu_0*i/2/pi*z/r/sqrt((a + r)**2 + z**2)*(-ellipk(k(r,z,a)) + (a**2 + r**2 + z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))

R = linspace(-10,10,100000)
Z = 0.01

t0 = time.time()
B1 = array([Br1(r,Z,0.5,10E3) for r in R])
t1 = time.time()
print(t1 - t0)

t0 = time.time()
B2 = array([Br2(r,Z,0.5,10E3) for r in R])
t1 = time.time()
print(t1 - t0)

plt.plot(R, B2 - B1)

Z = linspace(-10,10,10000)
R = 0.51

t0 = time.time()
B3 = array([Bz1(R,z,0.5,10E3) for z in Z])
t1 = time.time()
print(t1 - t0)

t0 = time.time()
B4 = array([Bz2(R,z,0.5,10E3) for z in Z])
t1 = time.time()
print(t1 - t0)

plt.plot(Z, B4 - B3)