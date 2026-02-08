import numpy as np
import matplotlib.pyplot as plt

num_k=400                       
kx=np.linspace(-0.5,0.5,num_k) 
ky=np.linspace(-0.5,0.5,num_k) 

G=3
G_vectors=[(n,m) for n in range(-G,G+1)
                     for m in range(-G,G+1)]
energies=[]

for n, m in G_vectors:
    E=(kx+n)**2+(kx+m)**2
    energies.append(E)

energies_sorted = np.sort(np.array(energies), axis=0)

num_bands=8

colors=['#1f77b4', '#7f7f7f', '#ff7f0e', '#2ca02c', '#FA8072',  '#003049', '#9467bd', '#e377c2']

plt.figure()
for i in range(num_bands):
    plt.plot(kx, energies_sorted[i], label=f'Band {i+1}', color=colors[i])

plt.axhline(y=7/(2*np.pi),color = "#2A9D8F",linestyle="--",label="Fermi energy")
plt.xlabel(r'$k$ (units of $2\pi/a$)')
plt.ylabel(r'$E$ (units of $\frac{\hbar^2}{2m}(\frac{2\pi}{a})^2$)')
plt.title('Empty lattice bands')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("6.1.png")
plt.show()
