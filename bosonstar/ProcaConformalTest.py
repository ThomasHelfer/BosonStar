from scipy.interpolate import interp1d
import os
import time

import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

def eqns(y, r, mu=1, cA4=0, GNewton=1):
    L, dLdr, Omega, a0, da0dr, a1 = y 
    
    dL2dr2 = dLdr**2/(2.*L) - (a1**4*cA4*GNewton*mu**2)/(4.*L) - (a1**2*GNewton*L*mu**2)/8. + (3*a0**4*cA4*GNewton*L**3*mu**2)/(4.*Omega**4) - (a1**2*GNewton*L)/(8.*Omega**2) + (a1*da0dr*GNewton*L)/(4.*Omega**2) - (da0dr**2*GNewton*L)/(8.*Omega**2) - (a0**2*a1**2*cA4*GNewton*L*mu**2)/(2.*Omega**2) - (a0**2*GNewton*L**3*mu**2)/(8.*Omega**2) - (2*dLdr)/r  

    dOdr = -((dLdr*Omega)/(L + dLdr*r)) - (a0**4*cA4*GNewton*L**3*mu**2*r)/(4.*Omega**3*(L + dLdr*r)) - (a1**2*GNewton*L*r)/(8.*Omega*(L + dLdr*r)) + (a1*da0dr*GNewton*L*r)/(4.*Omega*(L + dLdr*r)) - (da0dr**2*GNewton*L*r)/(8.*Omega*(L + dLdr*r)) - (a0**2*a1**2*cA4*GNewton*L*mu**2*r)/(2.*Omega*(L + dLdr*r)) + (a0**2*GNewton*L**3*mu**2*r)/(8.*Omega*(L + dLdr*r)) - (dLdr**2*Omega*r)/(2.*L*(L + dLdr*r)) + (3*a1**4*cA4*GNewton*mu**2*Omega*r)/(4.*L*(L + dLdr*r)) + (a1**2*GNewton*L*mu**2*Omega*r)/(8.*(L + dLdr*r))

    da1dr = (-4*a0*a1**2*cA4*L**2)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) + (8*a0*a1*cA4*da0dr*L**2)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (a0*L**4)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (a1*dLdr*L)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (da0dr*dLdr*L)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (4*a0**3*cA4*L**4)/(Omega**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (a1*dOdr*L**2)/(mu**2*Omega*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (da0dr*dOdr*L**2)/(mu**2*Omega*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (8*a1**3*cA4*dOdr*Omega)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (2*a1*dOdr*L**2*Omega)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) + (8*a1**3*cA4*dLdr*Omega**2)/(L*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (2*a1*L**2)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)*r) + (2*da0dr*L**2)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)*r)

    d2a0dr2 = da1dr + (a1*dLdr)/L - (da0dr*dLdr)/L + 4*a0*a1**2*cA4*mu**2 + a0*L**2*mu**2 - (4*a0**3*cA4*L**2*mu**2)/Omega**2 - (a1*dOdr)/Omega + (da0dr*dOdr)/Omega + (2*a1)/r - (2*da0dr)/r

    dydr = [dLdr, dL2dr2, dOdr, da0dr, d2a0dr2, da1dr]
    return dydr

#L, dLdr, Omega, a0, da0dr, a1 = y 
y0 = [1, 0, 0.8734, 0.165, 0, 0]
eps = 1e-10; R_max = 30; N = 10000
r = np.linspace(eps, R_max, N)
sol = spi.odeint(eqns, y0, r, atol = 1e-8, rtol = 1e-10)

plt.plot(r, sol[:, 0])
plt.title("L")
plt.figure()
plt.plot(r, sol[:, 2])
plt.title("Omega")
plt.figure()
plt.plot(r, sol[:, 3])
plt.title("a0")
plt.figure()
plt.plot(r, sol[:, 5])
plt.title("a1")
plt.show()