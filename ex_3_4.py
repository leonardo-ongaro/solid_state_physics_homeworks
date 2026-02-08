import numpy as np
import matplotlib.pyplot as plt

a = 1.0      
K = 1.0      
m = 1.0      

# Ex 3.b

def omega(kx, ky):
    return 2.0 * np.sqrt(K/m) * np.sqrt(
        np.sin(kx * a / 2.0)**2 + np.sin(ky * a / 2.0)**2
    )

Gamma = np.array([0.0, 0.0])
X     = np.array([np.pi/a, 0.0])
M     = np.array([np.pi/a, np.pi/a])

n_k = 200 

path = [
    np.linspace(Gamma, X, n_k),
    np.linspace(X, M, n_k),
    np.linspace(M, Gamma, n_k)
]

k_path = np.vstack(path)

omega_path = np.array([omega(k[0], k[1]) for k in k_path])

k_dist = np.zeros(len(k_path))
for i in range(1, len(k_path)):
    k_dist[i] = k_dist[i-1] + np.linalg.norm(k_path[i] - k_path[i-1])


plt.figure()
plt.plot(k_dist, omega_path, color='blue')

ts = [0,
         k_dist[n_k],
         k_dist[2*n_k],
         k_dist[-1]]

labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$']

for t in ts:
    plt.axvline(t, color='gray', linestyle='--', linewidth=0.8)

plt.xticks(ts, labels)
plt.ylabel(r'$\omega(k_x,k_y)$')
plt.xlabel('Path in ' + r'$\vec{k}$')
plt.title('Plot of ' + r'$\omega(k_x,k_y)$' + ' for a bidimensional square lattice')
plt.grid(True, which='both', axis='both', ls='--', linewidth=0.6, alpha=0.7)
plt.savefig("3.2.png")


#3.c
N = 2000 

kx_vals = np.linspace(-np.pi/a, np.pi/a, N)
ky_vals = np.linspace(-np.pi/a, np.pi/a, N)

kx_grid, ky_grid = np.meshgrid(kx_vals, ky_vals)


omega_vals = omega(kx_grid, ky_grid)

omega_flat = omega_vals.flatten()


n_bins = 100

hist, bin_edges = np.histogram(
    omega_flat,
    bins=n_bins,
    density=True  
)

centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure()
plt.plot(centers, hist, color='blue')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$D(\omega)$')
plt.title('Vibrational density of states')
plt.grid(True, ls='--', alpha=0.7)
plt.savefig("3.3.png")

# Ex 4.c
kB = 1.0
hbar = 1.0


domega = centers[1] - centers[0]

def specific_heat(T):
    x = centers / T
    integrand = D_omega = hist * (x**2) * np.exp(x) / (np.exp(x) - 1)**2
    return np.sum(integrand) * domega


T_vals = np.linspace(0.05, 6.0, 200)

C_vals = np.array([specific_heat(T) for T in T_vals])


plt.figure()
plt.plot(T_vals, C_vals, color='blue')
plt.xlabel(r'$T$')
plt.ylabel(r'$C_V$')
plt.title('Specific heat')
plt.grid(True, ls='--', alpha=0.7)
plt.savefig("4.5.png")
plt.show()




