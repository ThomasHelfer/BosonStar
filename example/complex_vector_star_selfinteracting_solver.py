from bosonstar.ComplexProcaSelfInteracting import Complex_Proca_Star
import numpy as np
# =====================
#  All important definitions
# =====================

# Physics definitions
f0 = 0.165          # central phi
D = 4.0             # Dimension (total not only spacial)
Lambda = 0.0        # Cosmological constant
mu = 1              # mass of the field
cA4 = 0.0           # selfinteraction of field
# Solver definitions
Rstart = 10.4
Rend = 50.00
deltaR = 2
N = 100000
# for G = 4 use sigma_guess = 0.779
GNewton = 1.0
sigma_guess = 0.8734

rtol = 1e-12
atol = 1e-14

verbose = 3
eps = 1e-10  # Small epsilon to avoid r =  0

# ====================================
#   Main routine
# ====================================

pewpew = Complex_Proca_Star(sigma_guess, f0, cA4, mu, GNewton, verbose, rtol,
        atol )
pewpew.print_parameters()

pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)

sol = pewpew.get_solution()

# =====================================
#   Output and plotting
# =====================================

pewpew.plot_solution()
pewpew.print_solution()
