'''
    Fig. 3
'''
import matplotlib.pyplot as plt
import numpy as np

dir = "/home_th/poley/Cascade_model_code_final/PythonFinal/Fig. 3, phase plots sig-gamma and mu-gamma planes/"

f = 1
fig, axs = plt.subplots(1, 2, figsize=(12*f, 5*f))

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


# general marker and other settings
lwidth = 3
msize = 800
malpha = 0.6
salpha = 1

# p = proportion of predator prey pairs
p = np.linspace(0, 1, 200)
gamma = np.cos(np.pi*p)

'''
    Panel (a)
'''

sig_rho1 = np.genfromtxt(dir + "sig-gamPlaneOppernu = 3.csv")
sig_rho1_plus = np.genfromtxt(dir + "sig-gamPlaneOppernu = 4.csv")
sig_nu1 = np.genfromtxt(dir + "sig-gamPlaneOpperrho = 3.csv")
sig_nu1_plus = np.genfromtxt(dir + "sig-gamPlaneOpperrho = 4.csv")

for idx in range(10):
    axs[0].arrow(p[20*idx+10], sig_rho1[20*idx+10]-0.08, 0, sig_rho1_plus[20*idx+10]-sig_rho1[20*idx+10], width = 0.008, color=color2)
    if idx < 6:
        axs[0].arrow(p[20*idx+10], sig_nu1[20*idx+10]+0.1, 0, sig_nu1_plus[20*idx+10]-sig_nu1[20*idx+10], width = 0.008, color=color3)
    if idx == 6:
        axs[0].arrow(p[20*idx+5], sig_nu1[20*idx+5]+0.15, 0, sig_nu1_plus[20*idx+5]-sig_nu1[20*idx+5], width = 0.008, color=color3)

sig = np.sqrt(2)/(1 + gamma)

axs[0].plot(p, sig_rho1, linewidth = lwidth, color = "black", linestyle="-", label=r"$\rho=4$")
axs[0].plot(p, sig_nu1, linewidth = lwidth, color = "black", linestyle="-", label=r"$\nu=4$")

axs[0].plot(p, sig, linewidth = lwidth, color = "black", label=r"$\rho=1, \nu=0$")

axs[0].set_xlabel(r"Proportion of PP Interactions")
axs[0].set_ylabel(r"$\sigma$")
axs[0].set_xticks([0.0, 1.0])
axs[0].set_yticks([0,3])
axs[0].set_ylim(0, 3)
axs[0].set_xlim(0, 1)


axs[0].fill_between(p, sig, color = color3)
axs[0].fill_between(p, sig, np.linspace(9, 9, 200), color = color2)
axs[0].fill_between(p[180:], np.linspace(9, 9, 20), color = color3)


axs[0].text(p[160]-0.16, sig_rho1[160]-1.1, r"Increasing $\rho$",bbox=dict(facecolor='white', alpha=0.5, boxstyle='round'))
axs[0].text(p[5]+0.02, sig_nu1_plus[5]+0.6, r"Increasing $\nu$",bbox=dict(facecolor='white', alpha=0.5, boxstyle='round'))

'''
    Panel (b)
'''

mu_p = np.genfromtxt(dir + "mu-gamplaneMinftynu=0rho=1.csv")
mu_p_nu3 = np.genfromtxt(dir + "mu-gamplaneMinftynu=3rho=1.csv")
mu_p_nu4 = np.genfromtxt(dir + "mu-gamplaneMinftynu=4rho=1.csv")
mu_p_rho3 = np.genfromtxt(dir + "mu-gamplaneMinftynu=0rho=3.csv")
mu_p_rho4 = np.genfromtxt(dir + "mu-gamplaneMinftynu=0rho=4.csv")

def mu_from_data(data, nu, rho):
    return (data[:, 2]*rho**2 + data[:, 1])/(data[:, 2]*rho**2 - data[:, 1])*nu

# setting nu = 0.05 and rho = 1.02 rather than nu = 0 and rho = 1
mu_p = mu_from_data(mu_p, 0.05, 1.02)
mu_p_nu3 = mu_from_data(mu_p_nu3, 3.00, 1.02)
mu_p_nu4 = mu_from_data(mu_p_nu4, 4.00, 1.02)
mu_p_rho3 = mu_from_data(mu_p_rho3, 0.05, 3.00)
mu_p_rho4 = mu_from_data(mu_p_rho4, 0.05, 4.00)

axs[1].plot(p[::-1], mu_p, linewidth = lwidth, color = "black", linestyle="-", label=r"$\rho=4$")
axs[1].plot(p[::-1], mu_p_nu3, linewidth = lwidth, color = "black", linestyle="-", label=r"$\rho=4$")
axs[1].plot(p[::-1], mu_p_rho3, linewidth = lwidth, color = "black", linestyle="-", label=r"$\rho=4$")

plt.fill_between(p[::-1], np.linspace(-9, -9, 200), mu_p, color = color3)
plt.fill_between(p[::-1], mu_p, np.linspace(9, 9, 200), color = color1, alpha = 0.85)



for idx in range(10):
    axs[1].arrow(p[::-1][20*idx+10], mu_p_nu3[20*idx+10]+0.16, 0, mu_p_nu4[20*idx+10]-mu_p_nu3[20*idx+10], width = 0.008, color=color3, head_length = 2*4.5*0.008)
    if idx < 8:
        axs[1].arrow(p[::-1][20*idx+10], mu_p_rho3[20*idx+10]-0.1, 0, mu_p_rho4[20*idx+10]-mu_p_rho3[20*idx+10], width = 0.008, color=color1, head_length = 2*4.5*0.008)
    if idx == 8:
        axs[1].arrow(p[::-1][20*idx+10], mu_p_rho3[20*idx+10]-0.1, 0, mu_p_rho4[20*idx+10]-mu_p_rho3[20*idx+10], width = 0.008, color=color3, head_length = 2*4.5*0.008)
    if idx == 9:
        axs[1].arrow(p[::-1][20*idx+10], mu_p_rho3[20*idx+10]+0.1, 0, mu_p_rho4[20*idx+10]-mu_p_rho3[20*idx+10], width = 0.008, color=color3, head_length = 2*4.5*0.008)


axs[1].set_xlabel(r"Proportion of PP Interactions")
axs[1].set_ylabel(r"$\mu$")
axs[1].set_xticks([0.0, 1.0])
axs[1].set_yticks([-3,0,3])
axs[1].set_ylim(-3, 3)
axs[1].set_xlim(0, 1)


axs[1].text(p[160]-0.15, mu_p_nu4[160]-1.6, r"Increasing $\rho$",bbox=dict(facecolor='white', alpha=0.5, boxstyle='round'))
axs[1].text(p[5]+0.02, mu_p_rho4[5]+1.7, r"Increasing $\nu$",bbox=dict(facecolor='white', alpha=0.5, boxstyle='round'))

axs[0].annotate(r'\textbf{(a)}', xy=(0.05, 0.05), xycoords='axes fraction')
axs[1].annotate(r'\textbf{(b)}', xy=(0.05, 0.05), xycoords='axes fraction')

# plt.savefig("M_to_infty.pdf")
plt.show()