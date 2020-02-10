from bosonstar.ComplexProcaStar import *
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib
matplotlib.use("agg")

# =====================
#  All imporntnat definitions
# =====================

# Physics defintions
f0 = 0.40         # centeral phi
D = 5.0             # Dimension (total not only spacial)
Lambda = -1.0      # Cosmological constant
mu = 1              # mass of the field
# Solver definitions
Rstart = 1.4
Rend = 100.00
deltaR = 1
N = 100000
sigma_guess = 0.32

verbose = 2
eps = 1e-10  # Small epsilon to avoid r \neq 0

# ====================================
#   Main routine
# ====================================

pewpew = Complex_Proca_Star(sigma_guess, f0, Lambda, mu, verbose)

pewpew.print_parameters()

pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)


# =====================================
#   Output and plotting
# =====================================


pewpew.plot_solution()
pewpew.print_solution()
