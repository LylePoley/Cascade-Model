'''
    Fig1. Cascade Matrix heatmap
'''

import numpy as np
import matplotlib.pyplot as plt
# from numba import jit
import matplotlib.colors as colors

# @jit(nopython = True)
def random_cascade_matrix(mu_L, mu_U, sigma_L, sigma_U, gamma, N = 250):
    J = np.zeros((N, N))

    for i in range(N):
        for j in range(i):
            z1 = np.random.normal(0, 1)
            z2 = np.random.normal(0, 1)

            J[i, j] = mu_L/N + sigma_L*z1/np.sqrt(N)
            J[j, i] = mu_U/N + sigma_U*(gamma*z1 + np.sqrt(1 - gamma**2)*z2)/np.sqrt(N)

    return J
fig, axs = plt.subplots(1, 3, figsize = (10, 3))
plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 12})
# for display
def block_matrix(mu, sizes):
    N = np.sum(sizes)
    c = np.hstack((np.array([0]), np.cumsum(sizes)[:-1]))

    J = np.zeros((N, N))

    for i in range(np.size(sizes)):
        for j in range(np.size(sizes)):
            for u in range(sizes[i]):
                for v in range(sizes[j]):
                    J[c[i] + u, c[j] + v] = mu[i, j]

    return J

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))

cmap = plt.get_cmap('viridis')
new_cmap = truncate_colormap(cmap, 0.12, 0.98)

mu1 = np.random.uniform(0, 1, 16).reshape((4, 4))
sizes1 = np.array([5, 4, 7, 3])

mu_cascade = np.array([[1, 2, 2, 2], [0, 1, 2, 2], [0, 0, 1, 2], [0, 0, 0, 1]])
sizes_cascade = np.array([7, 5, 6, 4])
J_cascade = block_matrix(mu_cascade, sizes_cascade)

J1 = block_matrix(mu1, sizes1)
axs[0].matshow(J1, cmap = "viridis")
axs[1].matshow(J_cascade, cmap = new_cmap)
axs[2].matshow(random_cascade_matrix(-1, 1, 0.0, 0.0, 0.0, 500), cmap = new_cmap)


# print(random_cascade_matrix)
# for i in range(3):
    # Turn off tick labels
axs[0].set_xticklabels([])
axs[0].set_xticks([])

axs[1].set_xticklabels([])
axs[1].set_xticks([])

axs[2].set_xticklabels([])
axs[2].set_xticks([])

axs[0].get_yaxis().set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[2].get_yaxis().set_visible(False)

axs[0].set_xlabel("Block Structured")
axs[1].set_xlabel("Cascade Model, Finite B")
axs[2].set_xlabel("Cascade Model, infinite B")
plt.tight_layout()

# panel labels
axs[0].annotate('(a)', xy=(0.05, 0.05), xycoords='axes fraction', zorder = 10, ma = 'center', bbox=dict(facecolor='white', alpha=0.6, boxstyle='round'))
axs[1].annotate('(b)', xy=(0.05, 0.05), xycoords='axes fraction', zorder = 10, ma = 'center', bbox=dict(facecolor='white', alpha=0.6, boxstyle='round'))
axs[2].annotate('(c)', xy=(0.05, 0.05), xycoords='axes fraction', zorder = 10, ma = 'center', bbox=dict(facecolor='white', alpha=0.6, boxstyle='round'))
fig.tight_layout()
plt.savefig("Interaction Matrix Cartoon.pdf")