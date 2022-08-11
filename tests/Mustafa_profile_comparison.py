import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

a1 = np.loadtxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_0.01/cA4_0.0/a1.dat")
r = np.loadtxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_0.01/cA4_0.0/rvals.dat")

def Mustafa_radial(r,c,a):
    return c * 0.76 * (a*r) /((1 + 0.0096 * (a*r) * (a*r))**16.0)

def Mustafa_radial2(r,k):
    a = np.sqrt(k) / (2*np.pi)
    c = -k**2
    return c * 0.76 * (a *r) /((1 + 0.0096 * (a*r) * (a*r))**16.0)

popt, pcov = curve_fit(Mustafa_radial, r, a1)
popt2, pcov2 = curve_fit(Mustafa_radial2, r, a1)

fit_radial = [Mustafa_radial(x,popt[0],popt[1]) for x in r] 
fit_radial2 = [Mustafa_radial2(x,popt2[0]) for x in r] 
plt.plot(r, a1)
plt.plot(r, fit_radial, "--")
plt.plot(r, fit_radial2, "r--")
plt.savefig("Mustafa_fit.png",dpi=600)
print("overal factor: ",popt[0])
print("r factor: ",popt[1])
print(popt[1]**2)