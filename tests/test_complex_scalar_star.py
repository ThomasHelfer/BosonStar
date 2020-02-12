import numpy as np 

from BosonStar.ComplexBosonStar import Complex_Boson_Star

# =====================
#  All imporntnat definitions
# =====================

def test():

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

    pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)
# =====================================
#   Output and plotting
# =====================================
    r, sol = pewpew.get_solution()
    pewpew.normalise_edelta()

    phi = sol[:, 2]
    m = sol[:, 1]
    e_pow_delta = 1 / sol[:, 0]

    e_pow_delta_ref = np.genfromtxt("tests/reference_scalar_star/Lambda_-0.2/Dim_5.0/phi0_0.4/edelta.dat")
    r_ref      = np.genfromtxt("tests/reference_scalar_star/Lambda_-0.2/Dim_5.0/phi0_0.4/rvals.dat")
    m_ref      = np.genfromtxt("tests/reference_scalar_star/Lambda_-0.2/Dim_5.0/phi0_0.4/m.dat")
    phi_ref    = np.genfromtxt("tests/reference_scalar_star/Lambda_-0.2/Dim_5.0/phi0_0.4/phi.dat")

    np.testing.assert_array_equal(m,m_ref)
    np.testing.assert_array_equal(e_pow_delta,e_pow_delta_ref)
    np.testing.assert_array_equal(phi,phi_ref)
    np.testing.assert_array_equal(r,r_ref)
