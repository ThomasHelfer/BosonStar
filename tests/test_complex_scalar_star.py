import numpy as np

from bosonstar.ComplexBosonStar import Complex_Boson_Star

# =====================
#  All important definitions
# =====================


def test():

    # Physics definitions
    phi0 = 0.40         # central phi
    D = 5.0             # Dimension (total not only spacial)
    Lambda = -0.2       # Cosmological constant
# Solver definitions
    Rstart = 3
    Rend = 50.00
    deltaR = 1
    N = 100000
    e_pow_minus_delta_guess = 0.4999
    verbose = 3
    eps = 1e-10  # Small epsilon to avoid r \neq 0

# ====================================
#   Main routine
# ====================================

    pewpew = Complex_Boson_Star(
        e_pow_minus_delta_guess, phi0, D, Lambda, verbose)

    pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)
# =====================================
#   Output and plotting
# =====================================
    sol = pewpew.get_solution()

    pewpew.normalise_edelta()

    # This presents the residue of Res = G_{\mu\nu} - T_{\mu\nu
    Einstein_residue = pewpew.check_Einstein_equation()

    threshold_MSE = 1e-6
    Nbuffer = 10
    Einstein_residue_MSE = \
        np.square(Einstein_residue[:, Nbuffer:-Nbuffer]).mean()
    assert(Einstein_residue_MSE < threshold_MSE)


test()
