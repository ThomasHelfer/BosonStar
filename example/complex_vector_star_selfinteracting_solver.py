from BosonStar.ComplexProcaSelfInteracting import Complex_Proca_Star

# =====================
#  All imporntnat definitions
# =====================

# Physics defintions
f0 = 0.165        # centeral phi
D = 4.0             # Dimension (total not only spacial)
Lambda = 0.0        # Cosmological constant
mu = 1              # mass of the field
cA4 = 0.5           # selfinteraction of field
# Solver definitions
Rstart = 10.4
Rend = 40.00
deltaR = 2
N = 100000
sigma_guess = 0.779

verbose = 2
eps = 1e-10  # Small epsilon to avoid r \neq 0

# ====================================
#   Main routine
# ====================================

pewpew = Complex_Proca_Star(sigma_guess, f0, cA4, mu, verbose)

pewpew.print_parameters()

pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)

# =====================================
#   Output and plotting
# =====================================

pewpew.plot_solution()
pewpew.print_solution()
