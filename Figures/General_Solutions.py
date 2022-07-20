import numpy as np
import scipy.special as sp
import scipy.integrate as integrate
import scipy.optimize as opt
import os

dir = os.getcwd()

# functions for the implementation of both simplififed cases of the theory, sigma_U = sigma_L and mu_U = mu_L

##############################################################################
################### W_k FUNCTIONS ############################################
##############################################################################

'''
    When x < -7 for each of w0, w1, w2, we replace the function by a 3rd or 4th
    order assymptotic expansion, this avoids numerical errors when considering
    ratios like w0(x)/w2(x).
'''


def _w(Delta):
    return 1/np.sqrt(2 * np.pi) * np.exp(-1/2 * Delta**2)


def w0(x):
    def _w0(_x):
        return 1/2*(1 + sp.erf(_x/np.sqrt(2)))

    return np.piecewise(x, [x < -7, x >= -7], [lambda x:-_w(x)/x*(1 - 1/x**2 + 3/x**4 - 15/x**6), _w0])


def w1(x):
    def _w1(_x):
        return np.exp(-_x**2/2)/np.sqrt(2*np.pi) + 1/2*_x*(1 + sp.erf(_x/np.sqrt(2)))

    return np.piecewise(x, [x < -7, x >= -7], [lambda x:_w(x)/x**2*(1 - 3/x**2 + 15/x**4), _w1])


def w2(x):
    def _w2(_x):
        return _x*np.exp(-_x**2/2)/np.sqrt(2*np.pi) + 1/2*(1 + _x**2)*(1 + sp.erf(_x/np.sqrt(2)))

    return np.piecewise(x, [x < -7, x >= -7], [lambda x:-2*_w(x)/x**3*(1 - 6/x**2 + 45/x**4), _w2])


def LogMean(x, y):
    def lm(x, y):
        return y*sp.exprel(sp.log1p(x/y - 1))

    return np.piecewise(x, [x == y, x != y], [lambda x: x, lambda x: lm(x, y)])

##############################################################################
################### CALCULATING [[X]] ########################################
##############################################################################

# the integral measure for all integrals


def dX(x, nu, sigma, rho, u):
    return u**2/(2*nu*u*w1(x) - 1/2*sigma**2*(rho**2 - 1/rho**2)*x*w2(x))

# condition for finding Delta star, the point at which the integral measure above diverges


def D_star(nu, sigma, rho, u):
    if(np.isclose(nu, 0.0)):
        return 0.0
    else:
        G = 1/4*sigma**2*(rho**2 - 1/rho**2)/(nu*u)
        if(G < 0.0):
            guess = -0.01

            def denom(x):
                return G - w1(x)/x*w2(x)

        if(G > 0.0):
            guess = 10.0

            def denom(x):
                return 1/G - x*w2(x)/w1(x)

        else:
            return np.inf

        return opt.fsolve(denom, [guess])[0]

# first coefficient in the laurent series of dX(x, nu, sigma, rho, u)


def a_(nu, sigma, rho, u, Ds=None):
    if(Ds == None):
        Ds = D_star(nu, sigma, rho, u)

    if(np.isinf(Ds)):
        return 0
    else:
        return 4*nu*u/((4*nu*u + 2*sigma**2*(rho**2 - 1/rho**2))*w0(Ds) - 3*sigma**2*(rho**2 - 1/rho**2)*w2(Ds))

# whats left over after the first laurent term has been taken out


def regularised_integrand(x, nu, sigma, rho, u, func, Ds=None):
    if(Ds == None):
        Ds = D_star(nu, sigma, rho, u)

    if(np.isinf(Ds)):
        return func(x)*dX(x, nu, sigma, rho, u)
    else:
        return func(x)*dX(x, nu, sigma, rho, u) - func(Ds)*a_(nu, sigma, rho, u, Ds)/(x - Ds)

# this should handle all the cases


def int_DDelta(nu, sigma, rho, u, D0, D1, func, Ds=None):
    if(Ds == None):
        Ds = D_star(nu, sigma, rho, u)

    # if the singularity is in the range of integration then cut it out, this is not an approximation
    if(np.minimum(D0, D1) <= Ds <= np.maximum(D0, D1)):
        def integrand_DD(x):
            return regularised_integrand(x, nu, sigma, rho, u, func, Ds)

        if(np.isinf(Ds)):
            return integrate.quad(integrand_DD, D0, D1)[0]
        else:
            return a_(nu, sigma, rho, u, Ds)*func(Ds)*np.log(np.abs((D1 - Ds)/(D0 - Ds))) + \
                integrate.quad(integrand_DD, D0, D1, limit=200,
                               points=[D0, D1, Ds])[0]

    # if not just integrate the function
    else:
        def integrand_DD(x):
            return func(x)*dX(x, nu, sigma, rho, u)

        return integrate.quad(integrand_DD, D0, D1)[0]

# the main integrals we're interested in, they correspond to [1], [w_0(\Delta)] and [w_2(\Delta)] up to a factor of u**2


def int_1(nu, sigma, rho, u, D0, D1, Ds=None):
    return int_DDelta(nu, sigma, rho, u, D0, D1, lambda x: 1, Ds)


def int_w0(nu, sigma, rho, u, D0, D1, Ds=None):
    return int_DDelta(nu, sigma, rho, u, D0, D1, w0, Ds)


def int_w2(nu, sigma, rho, u, D0, D1, Ds=None):
    return int_DDelta(nu, sigma, rho, u, D0, D1, w2, Ds)

######################################################################################
######################### GETTING u, \Delta_0, \Delta_1 ##############################
######################################################################################

# equations for obtaining u, Delta_0, Delta_1 from our model parameters


def cond_1(nu, sigma, rho, gamma, u, D0, D1, Ds=None):
    return 1 - int_1(nu, sigma, rho, u, D0, D1, Ds)


def cond_w0(nu, sigma, rho, gamma, u, D0, D1, Ds=None):
    return u*(1 - u) - gamma*sigma**2*int_w0(nu, sigma, rho, u, D0, D1, Ds)


def cond_w2(nu, sigma, rho, gamma, u, D0, D1, Ds=None):
    return u**2 - LogMean(rho**2, 1/rho**2)*sigma**2*int_w2(nu, sigma, rho, u, D0, D1, Ds)

# we have to write the system above like this so that scipy.optimize.fsolve() accepts it


def system_general(nu, sigma, rho, gamma, u, D0, D1):
    Ds = D_star(nu, sigma, rho, u)
    args = (nu, sigma, rho, gamma, u, D0, D1, Ds)
    return [cond_1(*args),
            cond_w0(*args),
            cond_w2(*args)]


def make_path(nu, sigma, rho, gamma, steps):

    return np.linspace(nu[0], nu[1], steps), \
        np.linspace(sigma[0], sigma[1], steps), \
        np.linspace(rho[0], rho[1], steps), \
        np.linspace(gamma[0], gamma[1], steps)


'''
    finding u, D0, D1 along some path through parameter space
    the seed to should be an approximate u, D0, D1 for path[0]
    seed is for nu = 0.1, sigma = 1, rho = 1.01, gamma = -0.5
    nu, sigma, rho and gamma should be initialised with make_path()
'''


def solve_along_path(nu, sigma, rho, gamma, seed=[1.31127639, 0.90283813, 0.90017523], printProgress=False, useAsSeed=False):
    if printProgress:
        print("u, D0, D1, [error in solution], D_star")

    rootData = np.zeros(shape=(nu.shape[0], 3))
    root = seed

    for i in range(nu.shape[0]):
        def bound_system(x):
            return system_general(nu=nu[i], sigma=sigma[i], rho=rho[i], gamma=gamma[i],
                                  u=x[0], D0=x[1], D1=x[2])

        root = opt.fsolve(bound_system, root)
        rootData[i] = root
        if printProgress:
            print(root, ' ', bound_system(root), ' ', D_star(
                nu[i], sigma[i], gamma[i], root[0]), ' ', i)

    if useAsSeed:
        return root
    else:
        return rootData[:, 0], rootData[:, 1], rootData[:, 2]

######################################################################################
######################### [[M]] and [[\phi]] #########################################
######################################################################################


def avgM(mu, nu, rho, D0, D1):
    inv = (D1*rho**2 + D0)/(D1*rho**2 - D0)*nu - mu

    return 1/inv


def avgPhi(sigma, gamma, u):
    return u*(1 - u)/(gamma*sigma**2)

######################################################################################
######################### ABUNDANCE DISTRIBUTIONS ####################################
######################################################################################

# helper function for abundance_distribution, equal to q_rho
def q_rho(nu, sigma, rho, u, D0, D, Ds=None):
    return 1/rho*np.exp(1/2*sigma**2*(rho**2 - 1/rho**2)*int_w2(nu, sigma, rho, u, D0, D, Ds))

# abundance distribution [[P(x)]]
def abundance_distribution(x, mu, nu, sigma, rho, u, D0, D1, Ds=None):
    # this prefactor is equal to [[M]]/[[w_1(\Delta)\sqrt{q_\rho}]]
    U_prefactor = u/2*((D1*rho + D0/rho) - (D1*rho - D0/rho)*mu/nu)

    def U(y):
        return U_prefactor/q_rho(nu, sigma, rho, u, D0, y, Ds)

    def integrand(y):
        return U(y)*np.exp(-1/2*(x*U(y) - y)**2)

    return u**2/np.sqrt(2*np.pi)*int_DDelta(nu, sigma, rho, u, D0, D1, integrand, Ds)

# get the function Delta(\alpha) once D0, D1 and u are known


def Delta(n, alpha, nu, sigma, rho, u, D0, D1):
    def func(x):
        return integrate.quad(n, 0, alpha)[0] - u**2*int_1(nu, sigma, rho, u, D0, x)

    return opt.fsolve(func, [np.minimum(D0, D1)])

# M(Delta)


def M_hierarchy(mu, nu, sigma, rho, u, D0, D1, D_alpha):
    prefactor = 1/2*u*((D1*rho+D0/rho) - (D1*rho - D0/rho)*mu/nu)

    return w1(D_alpha)*q_rho(nu, sigma, rho, u, D0, D_alpha)/prefactor


def phi_hierarchy(D_alpha):
    return w0(D_alpha)

######################################################################################
######################### LINEAR INSTABILITY #########################################
######################################################################################


def u_lin_instability(rho, gamma):
    return 1/(1 + gamma/LogMean(rho**2, 1/rho**2))


# seed tels you which of nu, sigma, rho, gamma is the seed value
# this function fills out the rest from this seed along the paramter path set out
'''
    finding a seed is hard, one strategy is to use the fact that we KNOW the seed when nu = 0,
    we know \Delta_0 = \Delta_1*rho**2, \Delta_0 \to \0, u = 1/(1 + \gamma/\ell) and u^2 = 1/2*\ell\sigma^2.

    The problem with this is that \Delta_0 \to 0, which makes the integrals diffiult numerically. 
    What seems to work better is just making some guess when \nu is large (say nu = 10), out here 
    the solution usually converges and you can then iterate inwards to \nu = 0, where the difficulties arise.

'''


def solve_along_lin_instability(nu=None, sigma=None, rho=None, gamma=None,
                                seed=None, D0_seed=None, D1_seed=None,
                                printProgress=True, steps=None, useAsSeed=False):

    if printProgress:
        print(seed + ", D0, D1, [error in solution], D_star")

    # rootData has three columns, the first is filled in with the values of the parameter we're looking for
    # the second two columns are for the corresponding D0 and D1 values
    rootData = np.zeros(shape=(steps, 3))

    # each of these if/elifs are basically the same, theyre adjusted for which parameter we are solving for and which are given
    # theres probably a better way of doing this but I don't know what it is
    if seed == "nu":

        # seed for the system
        root = [nu, D0_seed, D1_seed]

        for i in range(steps):
            # the system with u = u_lin_instability specifies that we are on the edge of linear instability
            def bound_system(x):
                return system_general(nu=x[0], sigma=sigma[i], rho=rho[i], gamma=gamma[i],
                                      u=u_lin_instability(rho[i], gamma[i]), D0=x[1], D1=x[2])

            # use the last loops solution as this loops guess and update
            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            # progress of the solution, problems usually occur when/if D_star is near either of D0 or D1.
            # D_star is the point at which our integral measure diverges, when it is not near D0 or D1
            # we can either ignore it or cut it out of our integrals
            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    root[0], sigma[i], rho[i], u_lin_instability(rho[i], gamma[i])), i)

    elif seed == "sigma":

        root = [sigma, D0_seed, D1_seed]

        for i in range(steps):
            def bound_system(x):
                return system_general(nu=nu[i], sigma=x[0], rho=rho[i], gamma=gamma[i],
                                      u=u_lin_instability(rho[i], gamma[i]), D0=x[1], D1=x[2])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], root[0], rho[i], u_lin_instability(rho[i], gamma[i])), i)

    elif seed == 'rho':

        root = [rho, D0_seed, D1_seed]

        for i in range(steps):
            def bound_system(x):
                return system_general(nu=nu[i], sigma=sigma[i], rho=x[0], gamma=gamma[i],
                                      u=u_lin_instability(x[0], gamma[i]), D0=x[1], D1=x[2])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], sigma[i], root[0], u_lin_instability(root[0], gamma[i])), i)

    elif seed == 'gamma':

        root = [gamma, D0_seed, D1_seed]

        for i in range(steps):
            def bound_system(x):
                return system_general(nu=nu[i], sigma=sigma[i], rho=rho[i], gamma=x[0],
                                      u=u_lin_instability(rho[i], x[0]), D0=x[1], D1=x[2])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], sigma[i], rho[i], u_lin_instability(rho[i], root[0])), ' ', i)

    if useAsSeed:
        return root
    else:
        return rootData[:, 0], rootData[:, 1], rootData[:, 2]

######################################################################################
######################### [[M]] -> INFINITY ##########################################
######################################################################################


# we introduce the extra condition with a new equation rather than substituting to avoid 0/0 errors
def Minfty_cond(mu, nu, rho, D0, D1):
    return D1*rho**2*(mu - nu) - D0*(mu + nu)


def system_Minfty(mu, nu, sigma, rho, gamma, u, D0, D1):
    original = system_general(nu, sigma, rho, gamma, u, D0, D1)
    return [original[0],
            original[1],
            original[2],
            Minfty_cond(mu, nu, rho, D0, D1)]


def solve_along_Minfty_path(mu, nu, sigma, rho, gamma, u, D0, D1, seed, steps, printProgress=True, useAsSeed=False):
    print(seed + ", u, D0, D1, [error in solution], D_star")

    rootData = np.zeros(shape=(steps, 4))
    if seed == "mu":
        root = [mu, u, D0, D1]
        for i in range(steps):

            def bound_system(x):
                return system_Minfty(mu=x[0], nu=nu[i], sigma=sigma[i], rho=rho[i], gamma=gamma[i],
                                     u=x[1], D0=x[2], D1=x[3])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], sigma[i], rho[i], root[1]), ' ', i)

    if seed == "nu":
        root = [nu, u, D0, D1]
        for i in range(steps):

            def bound_system(x):
                return system_Minfty(mu=mu[i], nu=x[0], sigma=sigma[i], rho=rho[i], gamma=gamma[i],
                                     u=x[1], D0=x[2], D1=x[3])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    root[0], sigma[i], rho[i], root[1]), ' ', i)

    if seed == "sigma":
        root = [sigma, u, D0, D1]
        for i in range(steps):

            def bound_system(x):
                return system_Minfty(mu=mu[i], nu=nu[i], sigma=x[0], rho=rho[i], gamma=gamma[i],
                                     u=x[1], D0=x[2], D1=x[3])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], root[0], rho[i], root[1]), ' ', i)

    if seed == "rho":
        root = [rho, u, D0, D1]
        for i in range(steps):

            def bound_system(x):
                return system_Minfty(mu=mu[i], nu=nu[i], sigma=sigma[i], rho=x[0], gamma=gamma[i],
                                     u=x[1], D0=x[2], D1=x[3])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], sigma[i], root[0], root[1]), ' ', i)

    if seed == "gamma":
        root = [gamma, u, D0, D1]
        for i in range(steps):

            def bound_system(x):
                return system_Minfty(mu=mu[i], nu=nu[i], sigma=sigma[i], rho=rho[i], gamma=x[0],
                                     u=x[1], D0=x[2], D1=x[3])

            root = opt.fsolve(bound_system, root)
            rootData[i] = root

            if printProgress:
                print(root, ' ', bound_system(root), ' ', D_star(
                    nu[i], sigma[i], rho[i], root[1]), ' ', i)

    if useAsSeed:
        return root
    else:
        return rootData[:, 0], rootData[:, 1], rootData[:, 2], rootData[:, 3]