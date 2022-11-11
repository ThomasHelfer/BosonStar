from cmath import e
from scipy.interpolate import interp1d
import os
import time

import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


class Complex_Proca_Star:

    sigma_guess = None
    _f0 = None
    _Dim = None
    _Lambda = None
    _mu = None
    _cA4 = None
    _GNewton = None

    # ------------------------------------------------------------
    # Scipy parameters
    # ------------------------------------------------------------
    _rtol = None
    _atol = None

    verbose = None
    path = None

    _omega = None
    _M = None
    _isrescaled = None
    _sigma_final = None
    __solution_array = None
    __solution_r_pos = None

    _finished_shooting = False

    def __init__(self, sigma_guess, f0, cA4=0.0, mu=1, GNewton=1, verbose=0,
            rtol = 1e-10, atol = 1e-10):

        self.sigma_guess = sigma_guess
        self._f0 = f0
        self._Dim = 4
        self._Lambda = 0
        self._mu = mu
        self._cA4 = cA4
        self._GNewton = GNewton
        self._atol = atol
        self._rtol = rtol

        # Will give more messages with increasing value
        self.verbose = verbose

        self.make_file()
        return None

    def print_parameters(self):
        print("----------------------------------------------------")
        print(r"The cosmological constant $\Lambda$ ", self._Lambda)
        print("The dimension of the problen         ", self._Dim)
        print(r"Central value of f_0                ", self._f0)
        print(r"Value of selfinteraction cA4        ", self._cA4)
        print(" Please cite https://arxiv.org/pdf/1805.09867.pdf ")
        print(" and future paper of Superradiant selfinteraction ")
        print("----------------------------------------------------")

    def eqns(self, y, r):
        """ Differential equation for scalar fields from arXiv:gr-qc/0309131

        Parameters:
            y (list with reals): current status vector ( a(r), alpha(r), f(r), pi(r) )
            r (real) : current position

        Returns:
            dydr (list with reals): derivative for y in r

        """
        Lambda = self._Lambda
        mu = self._mu
        L, dLdr, Omega, a0, da0dr, a1 = y 
        cA4 = self._cA4
        GNewton = self._GNewton
        # We defined pi = dfdx - g
        # Where sigma  = e^{-\delta}

        dL2dr2 = dLdr**2/(2.*L) - (a1**4*cA4*GNewton*mu**2)/(4.*L) - (a1**2*GNewton*L*mu**2)/8. + (3*a0**4*cA4*GNewton*L**3*mu**2)/(4.*Omega**4) - (a1**2*GNewton*L)/(8.*Omega**2) + (a1*da0dr*GNewton*L)/(4.*Omega**2) - (da0dr**2*GNewton*L)/(8.*Omega**2) - (a0**2*a1**2*cA4*GNewton*L*mu**2)/(2.*Omega**2) - (a0**2*GNewton*L**3*mu**2)/(8.*Omega**2) - (2*dLdr)/r  

        dOdr = -((dLdr*Omega)/(L + dLdr*r)) - (a0**4*cA4*GNewton*L**3*mu**2*r)/(4.*Omega**3*(L + dLdr*r)) - (a1**2*GNewton*L*r)/(8.*Omega*(L + dLdr*r)) + (a1*da0dr*GNewton*L*r)/(4.*Omega*(L + dLdr*r)) - (da0dr**2*GNewton*L*r)/(8.*Omega*(L + dLdr*r)) - (a0**2*a1**2*cA4*GNewton*L*mu**2*r)/(2.*Omega*(L + dLdr*r)) + (a0**2*GNewton*L**3*mu**2*r)/(8.*Omega*(L + dLdr*r)) - (dLdr**2*Omega*r)/(2.*L*(L + dLdr*r)) + (3*a1**4*cA4*GNewton*mu**2*Omega*r)/(4.*L*(L + dLdr*r)) + (a1**2*GNewton*L*mu**2*Omega*r)/(8.*(L + dLdr*r))

        da1dr = (-4*a0*a1**2*cA4*L**2)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) + (8*a0*a1*cA4*da0dr*L**2)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (a0*L**4)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (a1*dLdr*L)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (da0dr*dLdr*L)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (4*a0**3*cA4*L**4)/(Omega**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) + (a1*dOdr*L**2)/(mu**2*Omega*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (da0dr*dOdr*L**2)/(mu**2*Omega*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (8*a1**3*cA4*dOdr*Omega)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) - (2*a1*dOdr*L**2*Omega)/(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2) + (8*a1**3*cA4*dLdr*Omega**2)/(L*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)) - (2*a1*L**2)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)*r) + (2*da0dr*L**2)/(mu**2*(-4*a0**2*cA4*L**2 + (12*a1**2*cA4 + L**2)*Omega**2)*r)

        d2a0dr2 = da1dr + (a1*dLdr)/L - (da0dr*dLdr)/L + 4*a0*a1**2*cA4*mu**2 + a0*L**2*mu**2 - (4*a0**3*cA4*L**2*mu**2)/Omega**2 - (a1*dOdr)/Omega + (da0dr*dOdr)/Omega + (2*a1)/r - (2*da0dr)/r

        dydr = [dLdr, dL2dr2, dOdr, da0dr, d2a0dr2, da1dr]

        return dydr

    def shoot(self, Omega_at_zero, r, output=False):
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
        #L, dLdr, Omega, a0, da0dr, a1 = y 
        y0 = [1, 0, Omega_at_zero, self._f0, 0, 0]
        # Solve differential equaion
        sol = spi.odeint(self.eqns, y0, r, atol = self._atol, rtol = self._rtol)
        f_end = sol[-1, 3]

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
                print(
                    "Edelta at R = eps ",
                    sigma_guess_tmp[0],
                    " with Rmax ",
                    R_max)

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

    def normalise_OOmega(self, M):
        """ Extracts omega from OOmega. Then we perform the coordinate transformation t -> omega t
        """
        if self._M is None:
            print("----------------------------------------")
            print("WARNING: MASS EXTRACTION SHOULD HAPPEN BEFORE NORMALIZATION")
            print("----------------------------------------")
        if self._omega is None:
            if self.verbose >= 2:
                print("rescale r")
            r_end = self.__solution_r_pos[-1]
            end_value = (1-M/(2*r_end))/(1+M/(2*r_end))
            print("OOmega end value:", end_value)
            omega = end_value / self.__solution_array[-1, 2]
            self._omega = omega
            # Renormalising the sigma
            self.__solution_array[:, 2] *= omega
            self.__solution_array[:, 3] *= omega
            self.__solution_array[:, 4] *= omega
        else:
            print(" sigma has been already normalised ")

    def normalise_Lambda(self, factor):
        """ Rescale r so that Lambda->1 at infinity
        """
        if self._isrescaled is None:
            if self.verbose >= 2:
                print("Normalise sigma ")
            #Lambda_infty = self.__solution_array[-1,0]
            #print(Lambda_infty)
            self._isrescaled = True
            # Renormalising Lambda
            self.__solution_array[:, 0] *= 1/factor
            self.__solution_array[:, 1] *= 1/factor
            self.__solution_array[:, 4] *= 1/factor
            self.__solution_array[:, 5] *= 1/factor
            # redefining r
            self.__solution_r_pos *= factor

        else:
            print(" sigma has been already normalised ")

    def mass_extraction(self, fit_radius):
        """  
        Get the mass of the star from the fit of the metric; returns the mass 
        """
        if self._omega is not None or self._isrescaled is not None:
            print("----------------------------------------")
            print("WARNING: MASS EXTRACTION SHOULD HAPPEN BEFORE NORMALIZATION")
            print("----------------------------------------")

        edge_index = np.searchsorted(self.__solution_r_pos, fit_radius)
        Lambda_schwarz = self.__solution_array[:, 0][edge_index:]

        def fit_func(x,M,K):
            return K * (M/(2*x) + 1)**2
        def schw_OOmega(x,M):
            return (1-M/(2*x))/(1+M/(2*x))

        fit_params =opi.curve_fit(fit_func, self.__solution_r_pos[edge_index:], Lambda_schwarz)
        self._M = fit_params[0][0]*fit_params[0][1]


        self.normalise_OOmega(fit_params[0][0])
        self.normalise_Lambda(fit_params[0][1])
        OOmega_schwarz = self.__solution_array[:, 2][edge_index:]

        fit_Lambda = np.array([fit_func(x, fit_params[0][0]*fit_params[0][1], 1) for x in 
                        self.__solution_r_pos[edge_index:]])
        fit_OOmega = np.array([schw_OOmega(x, fit_params[0][0]*fit_params[0][1]) for x in 
                        self.__solution_r_pos[edge_index:]])

        _, (ax1, ax2) = plt.subplots(2, figsize = (10,10))
        ax1.plot(self.__solution_r_pos[edge_index:], OOmega_schwarz, "--", 
                    alpha=0.5, label = "cut solution")
        ax1.plot(self.__solution_r_pos, self.__solution_array[:, 2], "g",
                    alpha = 0.7, label = "full solution")
        ax1.plot(self.__solution_r_pos[edge_index:], fit_OOmega, label = "fit solution")
        ax1.legend()
        ax2.plot(self.__solution_r_pos[edge_index:], Lambda_schwarz, "--", 
                    alpha=0.5, label = "cut solution")
        ax2.plot(self.__solution_r_pos, self.__solution_array[:, 0], "g",
                    alpha = 0.7, label = "full solution")
        ax2.plot(self.__solution_r_pos[edge_index:], fit_Lambda, label = "fit solution")
        ax2.legend()
        path = self.get_path
        plt.savefig(self.path + "/fit_test.png")

        return fit_params


    def make_file(self):
        """ Creates Folder for current physics problem if they do not yet exist
        """
        name_Field = "vector_field_star_self_interacting"
        name_Lambda = "/Lambda_" + str(self._Lambda)
        name_Dim = "/Dim_" + str(self._Dim)
        name_cA4 = "/cA4_" + str(self._cA4)
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
        path += name_cA4
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
             solution_array (dictonary) : solution array for Rmax
        """
        if self.__solution_array is None or self.__solution_r_pos is None:
            print("----------------------------------------")
            print("WARNING: SHOOTING HAS NOT BEEN PERFORMED")
            print("----------------------------------------")
            return None
        else:
            if self._omega is None:
                self.normalise_OOmega()
            if self._isrescaled is None:
                self.normalise_Lambda()
            soldict = {
                "rpos": self.__solution_r_pos,
                "a1": self.__solution_array[:, 4],
                "da0dr": self.__solution_array[:, 3],
                "a0": self.__solution_array[:, 2],
                "m": self.__solution_array[:, 1],
                "sigma": self.__solution_array[:, 0],
                "omega": self._omega,
                "mass": self._M
            }
            return soldict

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
            if self._omega is None:
                self.normalise_OOmega()
            if self._isrescaled is None:
                self.normalise_Lambda()

            a1 = self.__solution_array[:, 5]
            da0dr = self.__solution_array[:, 4]
            a0 = self.__solution_array[:, 3]
            OOmega = self.__solution_array[:, 2]
            dLdr = self.__solution_array[:, 1]
            L = self.__solution_array[:, 0]
            r = self.__solution_r_pos
            omega = self._omega
            M = self._M

            D = self._Dim
            Lambda = self._Lambda
            mu = self._mu

            np.savetxt(self.path + "/omega.dat", [omega])
            np.savetxt(self.path + "/mass.dat", [M])
            np.savetxt(self.path + "/rvals.dat", r)
            np.savetxt(self.path + "/L.dat", L)
            np.savetxt(self.path + "/OOmega.dat", OOmega)
            np.savetxt(self.path + "/dLdr.dat", dLdr)
            np.savetxt(self.path + "/a0.dat", a0)
            np.savetxt(self.path + "/da0dr.dat", da0dr)
            np.savetxt(self.path + "/a1.dat", a1)

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
            if self._omega is None:
                self.normalise_OOmega()
            if self._isrescaled is None:
                self.normalise_Lambda()
            if self.verbose >= 1:
                print("Plotting started")
            if self.verbose >= 1:
                start = time.time()

            a1 = self.__solution_array[:, 5]
            da0dr = self.__solution_array[:, 4]
            a0 = self.__solution_array[:, 3]
            OOmega = self.__solution_array[:, 2]
            dLdr = self.__solution_array[:, 1]
            L = self.__solution_array[:, 0]
            r = self.__solution_r_pos
            omega = self._omega
            real_da0dr = np.array([ (a0[i+1]-a0[i])/(r[i+1]-r[i]) for i in range(len(a0)-1)])
            real_da0dr = np.append(real_da0dr,[0])

            def second_deriv(x,h):
                out = np.zeros(len(x))
                for i in range(len(x)-2):
                    out[i] = (x[i+2] -2*x[i+1] + x[i])/(h**2)
                return out
            psi = np.sqrt(np.abs(L))
            lap_tmp = second_deriv(r * psi, r[1]-r[0])
            lap_psi = lap_tmp / (r * r)

            # find 90 % radius of R
            Rguess = 0.01
            maxa0 = max(a0)
            a0_tmp_fun = interp1d(r, a0 - maxa0 * 0.1)
            root = opi.root(a0_tmp_fun, Rguess)
            R90 = root.x[0]

            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(10, 10))
            ax1.plot(r, L, 'b', )
            ax2.plot(r, OOmega, 'g')
            ax3.plot(r, a0, 'r')
            ax4.plot(r, a1, 'black')

            ax3.set_xlabel('t')

            ax1.set_ylabel(r'$\Lambda $')
            ax2.set_ylabel(r'$\Omega $')
            ax3.set_ylabel(r'$a0 (t)$')
            ax4.set_ylabel(r'$a1 (t)$')

            ax1.grid()
            ax2.grid()
            ax3.grid()

            plt.savefig(self.path + "/overview.png")

            plt.figure()
            plt.plot(real_da0dr, "--")
            plt.plot(da0dr, "r",alpha=0.5)
            plt.savefig(self.path + "/da0dr_consistency.png", dpi=600)

            if self.verbose >= 1:
                print("Plotting finished in ", time.time() - start, " sec")
