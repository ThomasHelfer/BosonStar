import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

amplitude = "0.05"
a1 = np.loadtxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_" + amplitude + "/cA4_0.0/a1.dat")
r = np.loadtxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_" + amplitude + "/cA4_0.0/rvals.dat")
omega = np.loadtxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_" + amplitude + "/cA4_0.0/omega.dat")

def Mustafa_radial(r,c,a):
    return c * 0.76 * (a*r) /((1 + 0.0096 * (a*r) * (a*r))**16.0)

def Mustafa_radial2(r,k):
    a = k
    c = -k**2 * np.sqrt(2)**1
    return c * 0.76 * (a *r) /((1 + 0.0096 * (a*r) * (a*r))**16.0)

popt, pcov = curve_fit(Mustafa_radial, r, a1)
popt2, pcov2 = curve_fit(Mustafa_radial2, r, a1)

fit_radial = [Mustafa_radial(x,popt[0],popt[1]) for x in r] 
fit_radial2 = [Mustafa_radial2(x,popt2[0]) for x in r] 
plt.plot(r, a1, "b",label="data")
#plt.plot(r, fit_radial, "--")
plt.plot(r, fit_radial2, "r--",label="fit")
#plt.plot(r, omega*a1,"b--")
#plt.show()
plt.legend()
plt.savefig("Mustafa_fit_f" + amplitude + ".png",dpi=600)
print("mu/m factor: ",popt2[0]**2)
np.savetxt("../vector_field_star_self_interacting/Lambda_0/Dim_4/f0_" + amplitude + "/cA4_0.0/factor.dat", popt2)
