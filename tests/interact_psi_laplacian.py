import matplotlib.pyplot as plt
import numpy as np
r = np.loadtxt("/home/zipeng/current_research/2022/Projects/ProcaStarHOC/Proca_isotropic/vector_field_star_self_interacting/Lambda_0/Dim_4/f0_0.165/cA4_0.0/rvals.dat")
Lambda = np.loadtxt("/home/zipeng/current_research/2022/Projects/ProcaStarHOC/Proca_isotropic/vector_field_star_self_interacting/Lambda_0/Dim_4/f0_0.165/cA4_0.0/L.dat")
def second_deriv(x,h):
    out = np.zeros(len(x))
    for i in range(len(x)-2):
        out[i] = (x[i+2] -2*x[i+1] + x[i])/(h**2)
    return out
psi = np.sqrt(np.abs(Lambda))
lap_tmp = second_deriv(r * psi, r[1]-r[0])
lap_psi = lap_tmp / (r * r)

plt.figure()
plt.plot(r[5:], lap_psi[5:])
plt.figure()
plt.plot(r, psi)
plt.show()