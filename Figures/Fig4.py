import General_Solutions as gs
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.optimize as opt

dir = "/home_th/poley/Cascade_model_code_final/PythonFinal/Fig. 4, Abundance distributions/"


# Plot settings
f = 2
fig, axs = plt.subplots(1, 3, figsize=(20*f, 5*f))

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 40})
N = 8

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
lwidth = 7
msize = 600
malpha = 0.8
salpha = 0.7

mu = 0.5
nu = 1.5
sigma = 0.15
rho = 1/1.5
gamma = -0.3

# the corresponding values of u, d0, d1 for the parameter choice above
# to see how to calculate this, see Examples.ipynb
u, d0, d1 = [ 1.00623867,  0.4621497,  24.10574192]

# data for the abundance distributions (for panels (a) and (b))
data_nu_withHierarchy = np.genfromtxt(dir + "AD_two_peaks.csv", skip_header = 7, names = None, delimiter = ',')
dataView_nu_withHierarchy = data_nu_withHierarchy.view(dtype = np.float64).reshape((-1, 500 + 9))

data_nu_withoutHierarchy = np.genfromtxt(dir + "AD_two_peaks_OB.csv", skip_header = 7, names = None, delimiter = ',')
dataView_nu_withoutHierarchy = data_nu_withoutHierarchy.view(dtype = np.float64).reshape((-1, 500 + 9))

'''
    Panel (a)
'''

m_withHierarchy = np.linspace(-5, 5, 500)
AD_withHierarchy = np.genfromtxt(dir + "AD_withHierarchy.csv")
axs[0].plot(m_withHierarchy, AD_withHierarchy, linewidth = lwidth, color = color1, zorder = 7)


m_withoutHierarchy = np.linspace(0, 4, 100)
AD_withoutHierarchy = np.genfromtxt(dir + "AD_withoutHierarchy.csv")
axs[0].plot(m_withoutHierarchy, AD_withoutHierarchy, linewidth = lwidth, color = color4, zorder = 5)

axs[0].hist(dataView_nu_withHierarchy[:, 8:].flatten(), bins = 50, density = True, alpha = malpha, color = color2, zorder = 6)
axs[0].hist(dataView_nu_withoutHierarchy[:, 8:].flatten(), bins = 50, density = True, alpha = malpha, color = color5, zorder = 4)


'''
    Panel (b)
'''

RANK_D_withHierarchy = 1-np.cumsum(AD_withHierarchy)[1:]*np.diff(m_withHierarchy)
RANK_D_withoutHierarchy = 1-np.cumsum(AD_withoutHierarchy)[1:]*np.diff(m_withoutHierarchy)

order_withHierarchy = np.sort(dataView_nu_withHierarchy[:, 9:], axis = 1)
order_withoutHierarchy = np.sort(dataView_nu_withoutHierarchy[:, 9:], axis = 1)

alpha = np.linspace(0, 1, 500)

axs[1].plot(RANK_D_withHierarchy, m_withHierarchy[1:], linewidth = lwidth, color = color1, zorder = 7)
axs[1].plot(RANK_D_withoutHierarchy, m_withoutHierarchy[1:], linewidth = lwidth, color = color4, zorder = 5)

axs[1].scatter(1-alpha[::30], np.average(order_withHierarchy, axis = 0)[::30], s = msize, alpha = salpha,
                    marker = 'D', color = color2, zorder = 4, label=r"$M$")
axs[1].scatter(1-alpha[::30], np.average(order_withoutHierarchy, axis = 0)[::30], s = msize, alpha = salpha, 
                    marker = 's', color = color5, zorder = 2, label=r"$M$")
'''
    Panel (c)
'''

a = np.linspace(0, 1)

D1 = gs.VDelta(lambda x : 1.0, a, nu, sigma, rho, u, d0, d1)
M_h1 = gs.M_hierarchy(mu, nu, sigma, rho, u, d0, d1, D1)

axs[2].plot(a, M_h1, linewidth = lwidth, color = color1, linestyle = '-', zorder = 4)
# 1.97353473 = value of M when rho = 1, nu = 0 and the other parameters are the same
axs[2].plot(a, np.linspace(1.97353473, 1.97353473, 50), linewidth = lwidth, color = color4, linestyle = '-', zorder = 3)


axs[2].scatter(alpha[::30], np.average(dataView_nu_withHierarchy[:, 9:], axis = 0)[::30], s = msize, alpha = salpha, 
                    marker = 'D', color = color2, zorder = 4, label=r"$M$")
axs[2].scatter(alpha[::30], np.average(dataView_nu_withoutHierarchy[:, 9:], axis = 0)[::30], s = msize, alpha = salpha, 
                    marker = 's', color = color5, zorder = 2, label=r"$M$")


'''
    axis titles, limits, etc
'''

axs[0].set_ylim(0, 1.5)
axs[0].set_xlim(0.1, 3.2)
axs[0].set_xticks([0, 1, 2, 3])
axs[0].set_yticks([0, 0.5, 1.0, 1.5])
axs[0].grid(True)
axs[0].set_xlabel(r"Abundance")
axs[0].set_ylabel(r"Frequency")


axs[1].set_xlim(-0.02, 1.02)
axs[1].grid(True)
axs[1].set_ylim(0.1, 3.0)
axs[1].set_yscale("log")
axs[1].set_xlabel(r"Rank")
axs[1].set_ylabel(r"Abundance")


axs[2].set_xlim(-0.02, 1.02)
axs[2].grid(True)
axs[2].set_xlabel(r"$r(\alpha)$")
axs[2].set_ylabel(r"Abundance $M(\alpha)$")
axs[2].set_ylim(0.1, 3.0)
axs[2].set_yscale("log")


axs[0].annotate(r'\textbf{(a)}', xy=(0.05, 0.05), xycoords='axes fraction', zorder = 10, bbox=dict(facecolor='white', alpha=1.0, boxstyle='round', pad = 0.15))
axs[1].annotate(r'\textbf{(b)}', xy=(0.05, 0.05), xycoords='axes fraction')
axs[2].annotate(r'\textbf{(c)}', xy=(0.05, 0.05), xycoords='axes fraction')

# plt.savefig("Abundance Distribution2.pdf")
