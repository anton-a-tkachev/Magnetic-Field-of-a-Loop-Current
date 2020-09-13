"""This module contains functions that calculate the magnetic field 
of a circular loop with electric current"""

from scipy.constants import pi, mu_0
from scipy.special import ellipk, ellipe
from numpy import sqrt

# Argument for the elliptic integrals
k = lambda r, z, a: 4*a*r/((a + r)**2 + z**2)

# Axial magnetic field
Hz = lambda r, z, a, i: i/2/pi/sqrt((a + r)**2 + z**2)*(ellipk(k(r,z,a)) + (a**2 - r**2 - z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))

# Radial magnetic field
Hr = lambda r, z, a, i: i/2/pi*z/r/sqrt((a + r)**2 + z**2)*(-ellipk(k(r,z,a)) + (a**2 + r**2 + z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a)))

# Axial magnetic flux
def Bz(r, z, a, i, mu = 1):
    return mu*mu_0*Hz(r, z, a, i)

# Radial magnetic flux
def Br(r, z, a, i, mu = 1):
    return mu*mu_0*Hr(r, z, a, i)
