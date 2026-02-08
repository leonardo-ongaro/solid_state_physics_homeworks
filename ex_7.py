import numpy as np
import matplotlib.pyplot as plt

a=1.0               
sigma=0.04*a      
N_G=2         
N_k=1000         

V0s=[0.5, 7.0] 
G=np.arange(-N_G, N_G + 1)  


def V_G(G_diff, V0, sigma, a):
    return -V0 * np.sqrt(np.pi) * (sigma / a) * np.exp(-(G_diff**2 * (np.pi * sigma / a)**2))

def delta(l,m):
    if l==m:
        return 1
    else:
        return 0

def A_matrix(k, V0):
    A=np.zeros((len(G),len(G)))
    for l in range(len(G)):
        for m in range(len(G)):
            A[l, m]=(k - G[m])**2 * delta(l,m)
            A[l, m]+=V_G(G[l]-G[m], V0, sigma, a)
    return A

k_v=np.linspace(-0.5, 0.5, N_k)  

for V0 in V0s:
    bands=[]
    for k in k_v:
        A=A_matrix(k, V0)
        E=np.linalg.eigh(A)[0]  
        bands.append(E)
    bands=np.array(bands)
    
  
    plt.figure(figsize=(6,4))
    for n in range(bands.shape[1]):
        plt.plot(k_v, bands[:, n], color='blue')
    plt.xlabel(r'$k$ (units of $2\pi/a$)')
    plt.ylabel(r'$E_n(k)$ (units of $\hbar^2/(2M) (2\pi/a)^2$)')
    plt.title(f'Band structure for V0 = {V0}')
    plt.grid(True)
    plt.savefig("7.1." + str(V0) + ".png")
    plt.show()
    
    E_f=bands.flatten()
    E_min, E_max=E_f.min(), E_f.max()
    E_vals=np.linspace(E_min, E_max, 1000)
    DOS, _ =np.histogram(E_f, bins=E_vals, density=True)
    
    plt.figure()
    plt.plot(E_vals[:-1], DOS, color='orange')
    plt.xlabel(r'$E$ (units of $\hbar^2/(2M) (2\pi/a)^2$)')
    plt.ylabel('DOS (a.u.)')
    plt.title(f'Density of States for V0 = {V0}')
    plt.grid(True)
    plt.savefig("7.2." + str(V0) + ".png")
    plt.show()
