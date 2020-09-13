"""This module contains functions that calculate the magnetic field 
of a circular loop with electric current"""

from scipy.constants import pi, mu_0
from scipy.special import ellipk, ellipe
from numpy import sqrt, array

# Argument for the elliptic integrals
k = lambda r,z,a: 4*a*r/((a + r)**2 + z**2)

# Axial magnetic flux
Bz = lambda r,z,a,i: mu_0*i/2/pi/sqrt((a + r)**2 + z**2)*(ellipk(k(r,z,a)) + (a**2 - r**2 - z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))

# Radial magnetic flux
Br = lambda r,z,a,i: mu_0*i/2/pi*z/r/sqrt((a + r)**2 + z**2)*(-ellipk(k(r,z,a)) + (a**2 + r**2 + z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))

B = lambda r,z,a,i: array([Br(r,z,a,i), Bz(r,z,a,i)])

# Magnetic field of a solenoid
def SBz(ID, OD, DZ, W, I, NR, NZ):
    pass