import os
import time

import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class Complex_Proca_Star:

    sigma_guess = None
    _f0 = None
    _Dim = None
    _Lambda = None
    _mu = None

    verbose = None
    path = None

    _omega = None
    _sigma_final = None
    __solution_array = None
    __solution_r_pos = None

    _finished_shooting = False

    def __init__(self, sigma_guess, f0, Lambda, mu=1, verbose=0):

        self.sigma_guess = sigma_guess
        self._f0 = f0
        self._Dim = 4
        self._Lambda = Lambda
        self._mu = mu

        # Will give more messages with increasing value
        self.verbose = verbose

        self.make_file()
        return None

    def print_parameters(self):
        print("----------------------------------------------------")
        print(r"The cosmological constant $\Lambda$ ", self._Lambda)
        print("The dimension of the problen         ", self._Dim)
        print(r"Central value of f_0                ", self._f0)
        print(" Please cite https://arxiv.org/pdf/1609.01735.pdf ")
        print("             https://arxiv.org/pdf/1508.05395.pdf ")
        print("----------------------------------------------------")

    def eqns(self, y, r):
        """ Differential equation for scalar fields from arXiv:gr-qc/0309131

        Parameters:
            y (list with reals): current status vector ( a(r), alpha(r), f(r), pi(r) )
            r (real) : current position

        Returns:
            dydr (list with reals): derivative for y in r

        """
        D = float(self._Dim)
        Lambda = self._Lambda
        mu = self._mu
        sigma, m, f, pi = y
        # We defined pi = dfdx - g
        # Where sigma  = e^{-\delta}

        F = (1 - 2 * m / r**(D - 3) - 2 * Lambda * r**2 / ((D - 2) * (D - 1)))

        g = -(pi / (F * mu**2 * sigma**2))

        dsigmadr = r * mu ** 2 * (sigma * g**2.0 + sigma**(-1) * f**2 / F**2)
        dmdr = r**(D - 2) * 0.5 * (pi**2 / sigma + mu **
                                   2 * (g**2 * F + sigma**(-2) * f**2 / F))

        dfdr = pi + g

        dpidr = (dsigmadr * F * pi * r - 2 * F * pi * sigma +
                 f * mu**2 * r * sigma) / (F * r * sigma)

        dydr = [dsigmadr, dmdr, dfdr, dpidr]

        return dydr

    def shoot(self, sigma_at_zero, r, output=False):
        """ Solves differential equation

        Parameters:
            sigma_at_zero (real): The lapse value guess at r = rmin
            r       (real array) : Radial points used for solver
            output  (bool)       : if True outputs whole solution array

        Returns:
            f_end (real):. The f value at r = rmax
            or
            sol     (real array) : array containg solution

        """

        # Define initial data vector
        y0 = [sigma_at_zero, 0, self._f0, 0]
        # Solve differential equaion
        sol = spi.odeint(self.eqns, y0, r)
        f_end = sol[-1, 2]

        if not output:
            return f_end
        else:
            return sol

    def radial_walker(self, r_start, r_end, delta_R, N, eps):
        """ Performs shooting for multiple radii rmax shooting process.

        Parameters:
            r_start (real) : first rmax for which shooting is performed
            r_end (real) : maximum rmax for which shooting is performed
            delta_R (real) : stelpsize
            N (real) : number of gridpoints

        Returns:
            alpha0 (real):. alpha0 for rmax
        """
        range_list = np.arange(r_start, r_end, delta_R)
        sigma_guess_tmp = self.sigma_guess

        if self.verbose >= 1:
            print("Shooting started")
        if self.verbose >= 1:
            start = time.time()

        for R_max in range_list:
            r = np.linspace(eps, R_max, N)

            def fun(x): return self.shoot(x, r)
            root = opi.root(fun, sigma_guess_tmp)
            sigma_guess_tmp = root.x

            if self.verbose >= 2:
                print("Edelta at R = eps ", sigma_guess_tmp[0], " with Rmax ", R_max)

        if self.verbose >= 1:
            print("Shooting finished in ", time.time() - start, "sec")

        self._finished_shooting = True
        output_solution = True
        self.__solution_r_pos = np.linspace(eps, r_end, N)
        self.__solution_array = self.shoot(
            sigma_guess_tmp[0],
            self.__solution_r_pos,
            output_solution)
        self._sigma_final = sigma_guess_tmp

        return sigma_guess_tmp[0]

    def normalise_sigma(self):
        """ Extractsomega for e_pow_delta by the coordinate transformation  t -> omega t

        Parameters:
            sol (real array) : were the sol[:,1] corresponds to edelta^(-1) and
                           and asymtotic value that does not go to 1
        Returns:
            omega (real): frequency of scalar field
            sol (real array) : sol array with fixed edelta
        """
        if self._omega is None:
            if self.verbose >= 2:
                print("Normalise sigma ")
            one_over_sigma = 1. / self.__solution_array[:, 0]
            N = len(one_over_sigma)
            omega = one_over_sigma[N - 1]
            one_over_sigma = one_over_sigma / omega
            self._omega = omega
            self.__solution_array[:, 0] = 1. / one_over_sigma
        else:
            print(" edelta has been already normalised ")

    def make_file(self):
        """ Creates Folder for current physics problem if they do not yet exist
        """
        name_Field = "vector_field_star"
        name_Lambda = "/Lambda_" + str(self._Lambda)
        name_Dim = "/Dim_" + str(self._Dim)
        name_Param = "/f0_" + str(self._f0)

        path = name_Field
        if not os.path.exists(path):
            os.mkdir(path)
        path += name_Lambda
        if not os.path.exists(path):
            os.mkdir(path)
        path += name_Dim
        if not os.path.exists(path):
            os.mkdir(path)
        path += name_Param
        if not os.path.exists(path):
            os.mkdir(path)
            if self.verbose >= 1:
                print("Create Folder with relative", path, ".")
        else:
            if self.verbose >= 1:
                print("Folder with path", path, "already exists.")

        self.path = path

    def get_path(self):
        """ return
              path (string): Realtive path used for outputs
        """
        if self.path is None:
            self.make_file()
        return self.path

    def get_solution(self):
        """return
             solution_array (real array) : solution array for Rmax
        """
        if self.__solution_array is None or self.__solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
            return None
        else:
            return self.__solution_r_pos, self.__solution_array

    def print_solution(self):
        """ Prints solution if shooting has been performed already

        """
        if self.path is None:
            self.make_file()
        if self.__solution_array is None or self.__solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
        else:
            if self.path is None:
                self.make_file()
            if self.verbose >= 2:
                print("Write out data")
            pi = self.__solution_array[:, 3]
            f = self.__solution_array[:, 2]
            m = self.__solution_array[:, 1]
            sigma = self.__solution_array[:, 0]
            r = self.__solution_r_pos
            if self._omega is None:
                self.normalise_sigma()
            omega = self._omega

            D = self._Dim
            Lambda = self._Lambda
            mu = self._mu

            F = (1 - 2 * m / r**(D - 3) - 2 *
                 Lambda * r**2 / ((D - 2) * (D - 1)))
            g = -(pi / (F * mu**2 * sigma**2))

            np.savetxt(self.path + "/omega.dat", [omega])
            np.savetxt(self.path + "/rvals.dat", r)
            np.savetxt(self.path + "/sigma.dat", sigma)
            np.savetxt(self.path + "/m.dat", m)
            np.savetxt(self.path + "/f.dat", f)
            np.savetxt(self.path + "/g.dat", g)

    def plot_solution(self):
        """ Prints solution if shooting has been performed already

        """
        if self.path is None:
            self.make_file()
        if self.__solution_array is None or self.__solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
        else:

            if self.verbose >= 1:
                print("Plotting started")
            if self.verbose >= 1:
                start = time.time()

            f = self.__solution_array[:, 2]
            m = self.__solution_array[:, 1]
            sigma = self.__solution_array[:, 0]
            r = self.__solution_r_pos

            # find 90 % radius of R
            Rguess = 0.01
            maxf = max(f)
            f_tmp_fun = interp1d(r, f - maxf * 0.1)
            root = opi.root(f_tmp_fun, Rguess)
            R90 = root.x[0]

            fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10, 10))
            ax1.plot(r, sigma, 'b', )
            ax2.plot(r, m, 'g')
            ax3.plot(r, f, 'r')

            ax3.set_xlabel('t')

            ax1.set_ylabel(r'$\sigma $')
            ax2.set_ylabel('$ m (t)$')
            ax3.set_ylabel(r'$f (t)$')

            ax1.set_xlim([0, R90 * 2])
            ax2.set_xlim([0, R90 * 2])
            ax3.set_xlim([0, R90 * 2])

            ax1.grid()
            ax2.grid()
            ax3.grid()

            plt.savefig(self.path + "/overview.png")

            if self.verbose >= 1:
                print("Plotting finished in ", time.time() - start, " sec")
