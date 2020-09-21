from magnetic import Ap, Bz
from scipy.constants import e, m_e, c
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

me = m_e
U = 30E3
g = 1 + e*U/me/c**2

AAA = lambda r,z: Ap(r, z - 0.200, 0.275, 50E3) + Ap(r, z - 0.8, 0.150, 10E3) + Ap(r, z - 2.500, 0.450, 96E3) + Ap(r, z - 3.500, 0.450, 480E3)

BBB = lambda r,z: Bz(r, z - 0.200, 0.275, 50E3) + Bz(r, z - 0.8, 0.150, 10E3) + Bz(r, z - 2.500, 0.450, 96E3) + Bz(r, z - 3.500, 0.450, 480E3)

P0 = e*AAA(0.020, 0.001)

Ueff = lambda r,z: P0**2/4/g/me/e*((2*e*AAA(r,z)/P0 - 1)**2 - 1)

zlist = np.linspace(0.001, 4, 1000)
rlist = np.linspace( 0.000001, 0.2, 1000)
R, Z = np.meshgrid(rlist, zlist)
F = Ueff(R,Z)
fig,ax=plt.subplots(1,1)
lvls = np.linspace(0, 1E3, 2)
cp = ax.contourf(Z, R, F, levels=list(lvls))
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Effective potential')
ax.set_xlabel('z (cm)')
ax.set_ylabel('r (cm)')
plt.show()

B1 = np.array([BBB(0, z) for z in zlist])

f = lambda r: r*BBB(r,3.5) + AAA(r,3.5) - P0/e
RL1 = fsolve(f,0.001)[0]
RL2 = np.sqrt(g**2 - 1)*me*c/e/np.max(B1)
print(RL1)
print(RL2)