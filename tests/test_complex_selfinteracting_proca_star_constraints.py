from bosonstar.ComplexProcaSelfInteracting import Complex_Proca_Star
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('agg')


def test():
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
    Rend = 40.00
    deltaR = 2
    N = 100000
    # for G = 4 use sigma_guess = 0.779
    GNewton = 1.0
    sigma_guess = 0.8734

    verbose = 3
    eps = 1e-10  # Small epsilon to avoid r =  0

    # ====================================
    #   Main routine
    # ====================================

    pewpew = Complex_Proca_Star(sigma_guess, f0, cA4, mu, GNewton, verbose)
    pewpew.print_parameters()

    pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)

    sol = pewpew.get_solution()

    r = sol['rpos']
    m = sol['m']
    a1 = sol['a1']
    a0 = sol['a0']
    da0dr = sol['da0dr']
    sigma = sol['sigma']
    omega = sol['omega']

    # Taking out values close to boundary, since numpy gradient doesn't know
    # about boundary conditions
    N_buffer = 10

    threshold_MSE = 1e-8

    dx = r[1] - r[0]
    d2a0dr2 = np.gradient(da0dr, dx)
    da1dr = np.gradient(a1, dx)
    dsigmadr = np.gradient(sigma, dx)
    dmdr = np.gradient(m, dx)

    # Calculating the EM Constraint equation
    constraint = np.sqrt(1 - (2 * m) / r) * (4 * a0**3 * cA4 * mu**2 * r**3 - a0 * mu**2 * (2 * m - r) * r * (-r + a1**2 * (8 * cA4 * m - 4 * cA4 * r)) * sigma**2 + (-2 * m + r)**2 * \
                         sigma * (-(da0dr * dsigmadr * r) + (2 * da0dr + d2a0dr2 * r - da1dr * omega * r) * sigma + a1 * (dsigmadr * omega * r - 2 * omega * sigma))) / (2. * r * (-2 * m + r)**2 * sigma**3)

    constraint_MSE = np.square(constraint[N_buffer:-N_buffer]).mean()

    # Calculating Hamiltonian constraint equation
    ham_constraint = -(-6 * a0**4 * cA4 * GNewton * mu**2 * r**4 + a0**2 * GNewton * mu**2 * (2 * m - r) * r**2 * (-r + a1**2 * (8 * cA4 * m - 4 * cA4 * r)) * sigma**2 + (-2 * m + r)**2 * sigma**2 * (da0dr**2 * GNewton * r**2 - 2 * a1 * da0dr * \
                       GNewton * omega * r**2 - 8 * dmdr * sigma**2 + 2 * a1**4 * cA4 * GNewton * mu**2 * (-2 * m + r)**2 * sigma**2 + a1**2 * GNewton * r * (-2 * m * mu**2 * sigma**2 + r * (omega**2 + mu**2 * sigma**2)))) / (4. * r**2 * (-2 * m + r)**2 * sigma**4)

    ham_constraint_MSE = np.square(ham_constraint[N_buffer:-N_buffer]).mean()

    assert(constraint_MSE < threshold_MSE)
    assert(ham_constraint_MSE < threshold_MSE)
