"""This module contains functions that calculate the magnetic field 
of a circular loop with electric current"""


from scipy.special import ellipk, ellipe
import numpy as np
import time

mu0 = 4*np.pi*1E-7

# Argument for the elliptic integrals
k = lambda r,z,a: 4*a*r/((a + r)**2 + z**2)

# Axial magnetic flux
Bz = lambda r,z,a,i: mu0*i/2/np.pi/np.sqrt((a + r)**2 + z**2)*(ellipk(k(r,z,a)) + (a**2 - r**2 - z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a))) if i != 0.0 else 0.0

# Radial magnetic flux
Br = lambda r,z,a,i: mu0*i/2/np.pi*z/r/np.sqrt((a + r)**2 + z**2)*(-ellipk(k(r,z,a)) + (a**2 + r**2 + z**2)/((a - r)**2 + z**2)*ellipe(k(r,z,a))) if r != 0.0 and i != 0 else 0.0

B = lambda r,z,a,i: np.array([Br(r,z,a,i), Bz(r,z,a,i)])

# Magnetic field of a solenoid
def BB(R, Z, ID, OD, DZ, W, I, NR, NZ):
    ZZ = np.linspace(-0.5*DZ, 0.5*DZ, NZ)
    RR = np.linspace(0.5*ID, 0.5*OD, NR)
    PP = np.array([[rr,zz] for rr in RR for zz in ZZ])
    BB = np.sum(np.array([B(R, Z - pp[1], pp[0], I*W/NR/NZ) for pp in PP]), axis=0)
    return BB

t0 = time.time()
B0 = BB(R = 0.1,     # r of the observation point (m)
        Z = 0.0,     # z of the observation point (m)
        ID = 0.250,  # winding inner diameter (m)
        OD = 0.300,  # winding outer diameter (m)
        DZ = 0.230,  # winding depth along z-axis (m)
        W = 400,     # number of turns
        I = 2.16E3,  # current in the coil (A)
        NR = 10,     # number of winding layers along r-axis
        NZ = 40)     # number of winding layers along z-axis
t1 = time.time()

print(B0)
print(t1 - t0)