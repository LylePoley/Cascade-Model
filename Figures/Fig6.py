'''
    Fig. 6 in the main text
'''
import matplotlib.pyplot as plt
import numpy as np
import General_Solutions as gs


dir = gs.dir + "/Fig. 2, phase plots rho-nu plane/"

f = 2
fig, axs = plt.subplots(1, 3, figsize=(18*f, 5.5*f))

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 40})
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
lwidth = 1
msize = 800
malpha = 0.6
salpha = 1


'''
    Panel (a)
'''

# middle plot, the one bucket stability plot in the sigma, mu plane
def _mu(gamma, D):
    return D*gs.w2(D)/(gs.w1(D)*(gs.w2(D) + gamma*gs.w0(D)))

def _sigma(gamma, D):
    return np.sqrt(gs.w2(D))/(gs.w2(D) + gamma*gs.w0(D))

D = np.concatenate((np.linspace(-1.0, 4.0, 200), np.linspace(4.01, 1000, 500)), axis = None)

# minfty line
minfty_mu = _mu(0.8, D)
minfty_sig = _sigma(0.8, D)
opper_sig = np.full(700,np.sqrt(2)/1.8)

# filling
axs[0].fill_between(np.linspace(-1.1, 0, 700), np.linspace(1, 1, 700), opper_sig, color = color2, alpha = salpha)
axs[0].fill_between(minfty_mu, np.minimum(minfty_sig, opper_sig), color = color3, alpha = salpha)
axs[0].fill_betweenx(minfty_sig, minfty_mu, np.full(700, 2), color = color1, alpha = salpha)
axs[0].fill_betweenx([0, 0.1], [1, 1], [2, 2], color = color1, alpha = salpha)

axs[0].set_xlim(-1.1, 1.1)
axs[0].set_ylim(0, 1.0)
axs[0].vlines(0, 0, 1.5)
axs[0].set_xlabel(r"$\mu$")
axs[0].set_ylabel(r"$\sigma$")

'''
    Panel (b)
'''

# bottom right, opper lines in the rho nu plane

opper_b = np.genfromtxt(dir + "oppersig=0.7.csv")
minfty_b = np.genfromtxt(dir + "minftysig=0.7mu=-0.5.csv")

# reverse to the axes can be aligned
opper_rho = opper_b[:, 0][::-1]
# interpolate so that the two curves have the same x coords, filling between two curves doesnt work without this
minfty_rho = np.concatenate((np.interp(np.linspace(-5, 1, 600), np.linspace(-5, 1, 500), minfty_b[:, 0]), np.full(400, 20)))

mask = (opper_rho <= minfty_rho)
axs[1].plot(np.linspace(-5, 5, 1000)[mask], opper_rho[mask], linewidth = lwidth, color = color2)
axs[1].plot(-np.linspace(-5, 5, 1000)[mask], 1/opper_rho[mask], linewidth = lwidth, color = color2)


axs[1].plot(np.linspace(-5, 5, 1000), minfty_rho, color = color1, linewidth = lwidth)
axs[1].plot(-np.linspace(-5, 5, 1000), 1/minfty_rho, color = color1, linewidth = lwidth)

# filling for different stability regions
# stable
axs[1].fill_between(np.linspace(-5, 5, 1000), np.minimum(minfty_rho, opper_rho), 
                                                1/np.minimum(minfty_rho, opper_rho)[::-1], color = color3, alpha = salpha)

# metastable
# # above
axs[1].fill_between(np.linspace(-5, 5, 1000), opper_rho, minfty_rho, where = minfty_rho >= opper_rho, color = color2, alpha = salpha)
# # below
axs[1].fill_between(-np.linspace(-5, 5, 1000), 1/opper_rho, 1/minfty_rho, where = minfty_rho >= opper_rho, color = color2, alpha = salpha)

# divergence
axs[1].fill_between(np.linspace(-5, 5, 1000), minfty_rho, np.full(1000, 20), color = color1, alpha = salpha)
axs[1].fill_between(-np.linspace(-5, 5, 1000), 1/minfty_rho, 1/np.full(1000, 20), color = color1, alpha = salpha)


# other settings 
axs[1].set_yscale("log")
axs[1].set_xlim(-5, 5)
axs[1].set_ylim(1/11, 11)
axs[1].hlines(1, -5, 5, zorder = 5)
axs[1].vlines(0, 0, 11, zorder = 5)

'''
    Panel (c)
'''

# Bottom left, Minfty lines in the rho nu plane
minfty_c = np.genfromtxt(dir + "minftysig=0.8mu=0.5.csv")
opper_c = np.genfromtxt(dir + "oppersig=0.8.csv")

minfty_nu = minfty_c[::-1, 0]
# interpolate so that the two curves have the same y coords, filling between two curves doesnt work without this
opper_nu = np.concatenate((np.interp(np.linspace(0.01, 5.156, 454), np.linspace(0.01, 5.155994326976409248e+00, 1200), opper_c[::-1, 0]), np.full(546, 10)))

axs[2].plot(minfty_nu, np.linspace(0.01, 11, 1000), color = color1, linewidth = lwidth)
axs[2].plot(-minfty_nu, 1/np.linspace(0.01, 11, 1000), color = color1, linewidth = lwidth)

# # filling for different stability regions
# stable
axs[2].fill_betweenx(np.linspace(0.01, 11, 1000), np.maximum(minfty_nu, opper_nu), np.full(1000, 6), alpha = salpha, color = color3)
axs[2].fill_betweenx(1/np.linspace(0.01, 11, 1000), np.minimum(-minfty_nu, -opper_nu), np.full(1000, -6), alpha = salpha, color = color3)

# metastable
axs[2].fill_betweenx(np.linspace(0.01, 11, 1000), minfty_nu, opper_nu, where=opper_nu >= minfty_nu, alpha = salpha, color = color2)
axs[2].fill_betweenx(1/np.linspace(0.01, 11, 1000), -minfty_nu, -opper_nu, where = opper_nu >= minfty_nu, alpha = salpha, color = color2)

# divergence
axs[2].fill_betweenx(np.linspace(0.01, 11, 1000), minfty_nu, 
                                                0, color = color1, alpha = salpha)
axs[2].fill_betweenx(1/np.linspace(0.01, 11, 1000), -minfty_nu, 
                                                0, color = color1, alpha = salpha)

axs[2].set_yscale("log")
axs[2].set_xlim(-5, 5)
axs[2].set_ylim(1/11, 11)
axs[2].hlines(1, -5, 5, zorder = 5)
axs[2].vlines(0, 0, 11, zorder = 5)

'''
    Other settings
'''
# points where the neighbouring plots are from
# [mu, sig] = [-0.5, 0.7], [0.5, 0.8]
axs[0].scatter([-0.5], [0.7], color = "yellow", s = 1000, edgecolor = 'black', linewidth = 5, zorder = 10)
axs[0].scatter([0.5], [0.8], color = "red", s = 1000, edgecolor = 'black', linewidth = 5, zorder = 10)
axs[1].scatter([0.0], [1.0], color = "yellow", s = 1000, edgecolor = 'black', linewidth = 5, zorder = 10)
axs[2].scatter([0.0], [1.0], color = "red", s = 1000, edgecolor = 'black', linewidth = 5, zorder = 10)

# labels and ticks
axs[0].set_xticks([-1, 0, 1])
axs[0].set_yticks([0.0, 0.5, 1.0])

axs[1].set_xticks([-5, 0, 5])
axs[1].set_yticks([0.1, 1, 10])
axs[1].set_xlabel(r"$\nu$")
axs[1].set_ylabel(r"$\rho$")

axs[2].set_xticks([-5, 0, 5])
axs[2].set_yticks([0.1, 1, 10])
axs[2].set_xlabel(r"$\nu$")
axs[2].set_ylabel(r"$\rho$")

# labels for different stability regions
axs[1].text(-4, 4.5, "Diverging\nAbundances", size=40,ma = 'center', bbox=dict(facecolor='white', alpha=0.65, boxstyle="round", pad = 0.25))
axs[1].text(0.6, 7.5, "Lin. Unstable", ma = 'center', bbox=dict(facecolor='white', alpha=0.65, boxstyle="round", pad = 0.25))
axs[1].text(2.5, 2, "Stable", ma = 'center', bbox=dict(facecolor='white', alpha=0.65, boxstyle="round", pad = 0.25))

# panel labels
axs[0].text(0.08, 0.06, '(a)', horizontalalignment='center', verticalalignment='center', transform = axs[0].transAxes)
axs[1].text(0.08, 0.06, '(b)', horizontalalignment='center', verticalalignment='center', transform = axs[1].transAxes)
axs[2].text(0.08, 0.06, '(c)', horizontalalignment='center', verticalalignment='center', transform = axs[2].transAxes)

fig.tight_layout()

# plt.savefig("Fig. 2: rho-nu plane phase plots.pdf")

