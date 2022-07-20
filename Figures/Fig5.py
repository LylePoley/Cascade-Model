import General_Solutions as gs
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate

dir = gs.dir + "/Fig. 5, survival distributions/"

def vDelta(alpha, nu, sigma, rho, u, d0, d1):
    rootData = np.zeros(np.size(alpha))
    root = d0

    def func(alpha, x):
        return integrate.quad(lambda x : 1, 0, alpha)[0] - gs.int_1(nu, sigma, rho, u, d0, x)

    for a in range(np.size(alpha)):
        root = opt.fsolve(lambda x : func(alpha[a], x), root)
        rootData[a] = root

    return rootData


# Plot settings
f = 1
fig, axs = plt.subplots(1, 1, figsize=(6*f, 5*f))

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 16})
N = 5

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0,1,N)))

color1 = next(plt.gca()._get_lines.prop_cycler)['color']
color2 = next(plt.gca()._get_lines.prop_cycler)['color']
color3 = next(plt.gca()._get_lines.prop_cycler)['color']
color4 = next(plt.gca()._get_lines.prop_cycler)['color']
color5 = next(plt.gca()._get_lines.prop_cycler)['color']
color6 = next(plt.gca()._get_lines.prop_cycler)['color']
color7 = next(plt.gca()._get_lines.prop_cycler)['color']
color8 = next(plt.gca()._get_lines.prop_cycler)['color']

# general marker and other settings
lwidth = 3.5
msize = 80
malpha = 0.8
salpha = 0.7

# start point
u, d0, d1 = [1.05134272, 3.32516367, 3.39267796]

# mu = -1.0
nu = 0.0
sigma = 0.3
rho = 0.99
gamma = -0.6
alpha = np.linspace(0, 1, 50)

sigma2 = 0.75
nu2 = 3.0
rho2 = 3.0

rho3 = 1/3.0

data = np.genfromtxt(dir + "SurvivalDistributions.csv", skip_header = 1)

u = data[:, 0]
d0 = data[:, 1]
d1 = data[:, 2]

d_ob = vDelta(alpha, nu, sigma2, rho, u[0], d0[0], d1[0])
d_rho = vDelta(alpha, nu, sigma2, rho3, u[1], d0[1], d1[1])
d_nu = vDelta(alpha, nu2, sigma2, rho, u[2], d0[2], d1[2])
d_both1 = vDelta(alpha, nu2, sigma2, rho3, u[3], d0[3], d1[3])
d_both2 = vDelta(alpha, nu2, sigma2, rho2, u[4], d0[4], d1[4])


axs.plot(alpha, gs.w0(d_ob), linewidth=lwidth, color = color1, zorder = 5)
axs.plot(alpha, gs.w0(d_rho), linewidth=lwidth, color = color2, zorder = 4)
axs.plot(alpha, gs.w0(d_nu), linewidth=lwidth, color = color3, zorder = 4)
axs.plot(alpha, gs.w0(d_both1), linewidth=lwidth, color = color4, zorder = 1)
axs.plot(alpha, gs.w0(d_both2), linewidth=lwidth, color = color5, zorder = 1)

axs.set_ylim(0, 1.02)
alpha = np.linspace(0, 1, 500)
for name, c, m in [("nu0rho1", color1, "d"), ("nu0rho3",color2, "s"), ("nu3rho1", color3, "o"), ("nu3rho0.33", color4, "v"), ("nu3rho3", color5, "^")]:
    data_nu = np.genfromtxt(dir + "nalphadistribution_" + name + ".csv", skip_header = 8, names = None, delimiter = ',')
    dataView_nu = data_nu.view(dtype = np.float64).reshape((-1, 500 + 9))
    if(name == "nu0rho3"):
        axs.scatter(1-alpha[::40], (np.count_nonzero(dataView_nu[:, 9:], axis = 0)/dataView_nu.shape[0])[::40], s = msize, alpha = salpha, marker = m, color = c, zorder = 2)
    else:
        axs.scatter(alpha[::40], (np.count_nonzero(dataView_nu[:, 9:], axis = 0)/dataView_nu.shape[0])[::40], s = msize, alpha = salpha, marker = m, color = c, zorder = 2)

axs.set_xlabel(r"$r(\alpha)$")
axs.set_ylabel(r"Survival rate $\phi(\alpha)$")
axs.set_xticks([0, 0.5, 1])
axs.set_yticks([0, 0.5, 1])
axs.grid(True)
axs.set_xlim(0, 1)

# plt.savefig("Fig. 5: Survival Probabilities.pdf")
