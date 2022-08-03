'''
    Fig S1. Proof that the equations for instability work
'''

import numpy as np
import matplotlib.pyplot as plt
import General_Solutions as gs

# working_directory = "/home_th/poley/Simulation_Infinite_Buckets_copy (5)/"
dir = gs.dir + "/Fig. S1, Data for Fig. 2 phase plots/"
fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# written so you just have to copy paste the files relative path
def getData(relativePath, skipHeader = 0, comment = "#", Names = True):
    return np.genfromtxt(dir + relativePath, delimiter = ',', names = Names, skip_header=skipHeader, comments=comment)


plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 20})
N = 4
# plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.CMRmap(np.linspace(0,1,N)))
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0,1,N)))

color1 = next(plt.gca()._get_lines.prop_cycler)['color']
color2 = next(plt.gca()._get_lines.prop_cycler)['color']
color3 = next(plt.gca()._get_lines.prop_cycler)['color']
color4 = next(plt.gca()._get_lines.prop_cycler)['color']
color5 = next(plt.gca()._get_lines.prop_cycler)['color']

# mu = 0.5 or -0.5
# sigma = 0.8 or 0.7
# gamma = 0.8

# general marker and other settings
lwidth = 5
msize = 800
malpha = 0.6
salpha = 1

for label, name, N, i in [[r"$\mu=-0.5, \sigma=0.5$","_sig=0.8_mu=0.5_gam=0.8", 200,1],[r"$\mu=-0.5, \sigma=0.7$","_sig=0.7_mu=-0.5_gam=0.8",200,0]]:


    data = getData("rho_nu_plane_stability"+name+".csv", 7, Names = None).reshape(-1, 5, N + 9)
    # print(data.shape)
    stabilityInfo = data[:, :, 5:9]
    stabilityInfo[:, :, -1] = 1-stabilityInfo[:, :, -1]

    minftyInfo = stabilityInfo[:, :, 0]
    opperInfo = stabilityInfo[:, :, 1:]

    opperInfo = np.prod(opperInfo, axis = 2)

    z0 = np.average(opperInfo, axis = 1)
    z1 = np.average(minftyInfo, axis = 1)

    Z = 0.3*z0 + 0.1
    Z[z1 <= 0.9] = 0.0

    paramInfo = data[:, :, :6]
    paramInfo = np.average(paramInfo, axis = 1)

    nu = paramInfo[:, 1]
    rho = paramInfo[:, 3]

    axs[i].scatter(nu, rho, c = Z, s = 25)


# bottom right, opper lines in the rho nu plane

opper_b = np.genfromtxt(dir + "oppersig=0.7.csv")
minfty_b = np.genfromtxt(dir + "minftysig=0.7mu=-0.5.csv")

# reverse to the axes can be aligned
opper_rho = opper_b[:, 0][::-1]
# interpolate so that the two curves have the same x coords, filling between two curves doesnt work without this
minfty_rho = np.concatenate((np.interp(np.linspace(-5, 1, 600), np.linspace(-5, 1, 500), minfty_b[:, 0]), np.full(400, 20)))

mask = (opper_rho <= minfty_rho)
axs[0].plot(np.linspace(-5, 5, 1000)[mask], opper_rho[mask], linewidth = lwidth, color = color2)
axs[0].plot(-np.linspace(-5, 5, 1000)[mask], 1/opper_rho[mask], linewidth = lwidth, color = color2)


axs[0].plot(np.linspace(-5, 5, 1000), minfty_rho, color = color1, linewidth = lwidth)
axs[0].plot(-np.linspace(-5, 5, 1000), 1/minfty_rho, color = color1, linewidth = lwidth)

# other settings 
axs[0].set_yscale("log")
axs[0].set_xlim(-5, 5)
axs[0].set_ylim(1/11, 11)
axs[0].hlines(1, -5, 5, zorder = 5)
axs[0].vlines(0, 0, 11, zorder = 5)

'''
    Panel (c)
'''

# Bottom left, Minfty lines in the rho nu plane
minfty_c = np.genfromtxt(dir + "minftysig=0.8mu=0.5.csv")
opper_c = np.genfromtxt(dir + "oppersig=0.8.csv")

minfty_nu = minfty_c[::-1, 0]
# interpolate so that the two curves have the same y coords, filling between two curves doesnt work without this
opper_nu = np.concatenate((np.interp(np.linspace(0.01, 5.156, 454), np.linspace(0.01, 5.155994326976409248e+00, 1200), opper_c[::-1, 0]), np.full(546, 10)))

axs[1].plot(minfty_nu, np.linspace(0.01, 11, 1000), color = color1, linewidth = lwidth)
axs[1].plot(-minfty_nu, 1/np.linspace(0.01, 11, 1000), color = color1, linewidth = lwidth)


axs[1].set_yscale("log")
axs[1].set_xlim(-5, 5)
axs[1].set_ylim(1/11, 11)
axs[1].hlines(1, -5, 5, zorder = 5)
axs[1].vlines(0, 0, 11, zorder = 5)

'''
    Other settings
'''

# labels and ticks
axs[0].set_xticks([-5, 0, 5])
axs[0].set_yticks([0.1, 1, 10])
axs[0].set_xlabel(r"$\nu$")
axs[0].set_ylabel(r"$\rho$")

axs[1].set_xticks([-5, 0, 5])
axs[1].set_yticks([0.1, 1, 10])
axs[1].set_xlabel(r"$\nu$")
axs[1].set_ylabel(r"$\rho$")


# panel labels
axs[0].text(0.08, 0.06, '(a)', horizontalalignment='center', verticalalignment='center', transform = axs[0].transAxes)
axs[1].text(0.08, 0.06, '(b)', horizontalalignment='center', verticalalignment='center', transform = axs[1].transAxes)

fig.tight_layout()

# plt.savefig("Phase Diagram Proof.pdf")