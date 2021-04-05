import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

path_to_script = "../example/complex_vector_star_selfinteracting_solver.py"

# executing the script
exec(open(path_to_script).read())

print(sol.keys())

r = sol['rpos']
m = sol['m']
a1 = sol['a1']
a0 = sol['a0']
da0dr = sol['da0dr']
sigma = sol['sigma']
omega = sol['omega']

dx = r[1]-r[0]
da0dr = np.gradient(a0, dx)
d2a0dr2 = np.gradient(da0dr, dx)
da1dr = np.gradient(a1, dx)
dsigmadr = np.gradient(sigma, dx)
dmdr =  np.gradient(m, dx)
GNewton = 0

constraint = np.sqrt(1 - (2*m)/r)*(4*a0**3*cA4*mu**2*r**3 - a0*mu**2*(2*m - r)*r*(-r + a1**2*(8*cA4*m - 4*cA4*r))*sigma**2 + (-2*m + r)**2*sigma*(-(da0dr*dsigmadr*r) + (2*da0dr + d2a0dr2*r - da1dr*omega*r)*sigma + a1*(dsigmadr*omega*r - 2*omega*sigma)))/(2.*r*(-2*m + r)**2*sigma**3)

hamconstraint = -(-6*a0**4*cA4*GNewton*mu**2*r**4 + a0**2*GNewton*mu**2*(2*m - r)*r**2*(-r + a1**2*(8*cA4*m - 4*cA4*r))*sigma**2 + (-2*m + r)**2*sigma**2*(da0dr**2*GNewton*r**2 - 2*a1*da0dr*GNewton*omega*r**2 - 8*dmdr*sigma**2 + 2*a1**4*cA4*GNewton*mu**2*(-2*m + r)**2*sigma**2 + a1**2*GNewton*r*(-2*m*mu**2*sigma**2 + r*(omega**2 + mu**2*sigma**2))))/(4.*r**2*(-2*m + r)**2*sigma**4)

plt.close()
plt.plot(r[100:],np.abs(hamconstraint[100:]))
plt.xlim([0,12])
#plt.yscale("log")
plt.savefig("test.png")
