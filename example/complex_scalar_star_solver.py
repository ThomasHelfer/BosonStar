import os

from BosonStar.ComplexBosonStar import Complex_Boson_Star

# =====================
#  All imporntnat definitions
# =====================

# Physics defintions
phi0 = 0.40         # centeral phi
D = 5.0             # Dimension (total not only spacial)
Lambda = -0.2       # Cosmological constant
# Solver definitions
Rstart = 3
Rend = 50.00
deltaR = 1
N = 100000
e_pow_minus_delta_guess = 0.4999
verbose = 1
eps = 1e-10  # Small epsilon to avoid r \neq 0

# ====================================
#   Main routine
# ====================================

pewpew = Complex_Boson_Star(e_pow_minus_delta_guess, phi0, D, Lambda, verbose)

pewpew.print_parameters()

alpha0 = pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)


# =====================================
#   Output and plotting
# =====================================
r, sol = pewpew.get_solution()

pewpew.normalise_edelta()

# ===============================
path = pewpew.get_path()
pewpew.plot_solution()
pewpew.print_solution()
