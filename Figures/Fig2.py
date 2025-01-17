'''
    Fig 2. in the main text
'''

# general marker and other settings
import numpy as np
import matplotlib.pyplot as plt
import General_Solutions as gs

# set this to whatever folder the Fig.1 data is in
dir = gs.dir + "/Fig. 2, M and phi against nu data/"

lwidth = 8
msize = 800
malpha = 0.6
pathlen = 1000

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 40})
fig, axs = plt.subplots(1, 2, figsize=(32,12))


nu = 5.0
sigma = 0.3
rho = 1.5
mu = 0.958

gamma1 = -1
gamma2 = 1
gamma3 = 0

'''
    theory for left plot
'''
'''
    paths1, 2, 3, 4 each contain three columns, the first
    is G = sigma**2 (rho**2 - 1/rho**2)/(4 * nu * u)
    the second and third are \Delta_0 and \Delta_1 respectively.
'''
theoryData1 = np.genfromtxt(dir + "path3.csv") # gamma = -1
theoryData2 = np.genfromtxt(dir + "path4.csv")
theoryData3 = np.genfromtxt(dir + "path2.csv")

nu1 = gs.make_path([0.01, nu], [sigma, sigma], [rho, rho], [gamma1, gamma1], pathlen+500)[0]
nu2 = gs.make_path([0.01, nu], [sigma, sigma], [rho, rho], [gamma2, gamma2], pathlen+200)[0]
nu3 = gs.make_path([0.01, nu], [sigma, sigma], [rho, rho], [gamma3, gamma3], pathlen+600)[0]

M1 = gs.avgM(mu, nu1[::4], rho, theoryData1[::4, 1], theoryData1[::4, 2])
skip = 100
M2 = gs.avgM(mu, nu2[skip:], rho, theoryData2[skip:, 1], theoryData2[skip:, 2])
M3 = gs.avgM(mu, nu3[:], rho, theoryData3[:, 1], theoryData3[:, 2])

axs[0].plot(nu1[::4], M1, c = "darkblue", linewidth = lwidth)
axs[0].plot(nu2[skip:], M2, c = "darkred", linewidth = lwidth)
axs[0].vlines(0.4, 0, 100, colors = 'darkred', alpha = 0.8, linestyles = 'dashed', linewidth = 1.5)
axs[0].plot(nu3[:], M3, c = "darkgreen", linewidth = lwidth)

'''
    data for left plot
'''

simData1 = np.genfromtxt(dir + "M_to_infinity_gam-1.csv", names = None, skip_header = 7, comments = "#", delimiter = ",").reshape(-1, 40, 259)
simData2 = np.genfromtxt(dir + "M_to_infinity.csv", names = None, skip_header = 7, comments = "#", delimiter = ",").reshape(-1, 40, 259)
simData3 = np.genfromtxt(dir + "M_to_infinity_gam-0.csv", names = None, skip_header = 7, comments = "#", delimiter = ",").reshape(-1, 40, 259)

avgSimData1 = np.average(simData1, axis = 1)
avgSimData2 = np.average(simData2, axis = 1)
avgSimData3 = np.average(simData3, axis = 1)

nu1_data = avgSimData1[:, 1]
nu2_data = avgSimData2[:, 1]
nu3_data = avgSimData3[:, 1]

mu1_data = avgSimData1[:, 0]
mu2_data = avgSimData2[:, 0]
mu3_data = avgSimData3[:, 0]

M1 = np.average(avgSimData1[:, 9:], axis = 1)
M2 = np.average(avgSimData2[:, 9:], axis = 1)
M3 = np.average(avgSimData3[:, 9:], axis = 1)


axs[0].scatter(nu1_data, M1, marker = "s", c = "darkblue", s = msize, zorder=3, alpha = malpha, label = r"$\gamma = -1$")
axs[0].scatter(nu2_data, M2, marker = "D", c="darkred", s = msize, zorder=3, alpha = malpha, label = r"$\gamma = 1$")
axs[0].scatter(nu3_data, M3, marker = "o", c = "darkgreen", s = msize, zorder=3, alpha = malpha, label = r"$\gamma = 0$")



'''
    theory for right plot
'''

# small bits of jibbery, avoiding using the limiting form for the survival rate when gamma = 0 by using gamma = 0.01
# the real deal is indistinguishable
gamma3 += 0.01
theoryData3 = np.genfromtxt(dir + "path1.csv")
nu3 = gs.make_path([0.01, nu], [sigma, sigma], [rho, rho], [gamma3, gamma3], pathlen + 400)[0]

G1 = theoryData1[:, 0]
u1 = sigma**2*(rho**2 - 1/rho**2)/(4*nu1*G1)
phi1 = u1*(1-u1)/(gamma1*sigma**2)

G2 = theoryData2[:, 0]
u2 = sigma**2*(rho**2 - 1/rho**2)/(4*nu2*G2)
phi2 = u2*(1-u2)/(gamma2*sigma**2)

G3 = theoryData3[:, 0]
u3 = sigma**2*(rho**2 - 1/rho**2)/(4*nu3*G3)
phi3 = u3*(1-u3)/(gamma3*sigma**2)

axs[1].plot(nu1[(phi1>0.995) | (nu1 > 0.5)], phi1[(phi1>0.995) | (nu1 > 0.5)], label = r"$\gamma = -1$", c = "darkblue", linewidth = lwidth) 
axs[1].plot(nu2[nu2>0.4], phi2[nu2>0.4], label = r"$\gamma = 1$", c = "darkred", linewidth = lwidth) 
axs[1].plot(nu3, phi3, label = r"$\gamma = 0$", c = "darkgreen", linewidth = lwidth)

axs[1].vlines(0.4, 0, 0.996, colors = 'darkred', alpha = 0.8, linestyles = 'dashed', linewidth = 1.5)

'''
    simulated data for right plot
'''

phi1 = np.count_nonzero(simData1[:, :, 9:], axis = (1, 2))/(np.count_nonzero(simData1[:, :, 5] != 0, axis = 1)*250)
phi2 = np.count_nonzero(simData2[:, :, 9:], axis = (1, 2))/(np.count_nonzero(simData2[:, :, 5] != 0, axis = 1)*250)
phi3 = np.count_nonzero(simData3[:, :, 9:], axis = (1, 2))/(np.count_nonzero(simData3[:, :, 5] != 0, axis = 1)*250)

nu1_data = np.average(simData1, axis = 1)[:, 1]
nu2_data = np.average(simData2, axis = 1)[:, 1]
nu3_data = np.average(simData3, axis = 1)[:, 1]

axs[1].scatter(nu1_data, phi1, c = "darkblue", alpha = malpha, s = msize, marker = "s", zorder = 3)
axs[1].scatter(nu2_data, phi2, c = "darkred", alpha = malpha, marker = "D", s = msize, zorder = 3)
axs[1].scatter(nu3_data, phi3, c = "darkgreen", alpha = malpha, s = msize, zorder = 3)

'''
    options for the first plot
'''

axs[0].legend(prop={'size': 32})
axs[0].set_xlim(0, 1.05)
axs[0].set_yscale("log")
axs[0].set_ylim(1, 100.0)
axs[0].grid("true")

axs[0].set_ylabel(r'Average Abundance $[\![ M ]\!]$')
axs[0].set_xlabel(r'$\nu$')

'''
    plot settings for right plot
'''

axs[1].set_xlim(0, 1.05)
axs[1].set_ylim(0.979, 1)
axs[1].set_yticks([0.98, 0.99, 1.00])
axs[1].grid("true")


axs[1].set_ylabel(r'Fraction of Survivors $[\![ \phi ]\!]$')
axs[1].set_xlabel(r'$\nu$')

axs[0].text(0.05, 0.08, '(a)', horizontalalignment='center', verticalalignment='center', transform = axs[0].transAxes)
axs[1].text(0.05, 0.08, '(b)', horizontalalignment='center', verticalalignment='center', transform = axs[1].transAxes)
fig.tight_layout()
# plt.savefig("Fig. 2: theory working demonstration.pdf")