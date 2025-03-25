import numpy as np
import matplotlib.pyplot as plt

# Parâmetros da malha (ajuste conforme seu caso)
nx = 20
ny = 10
npx = nx + 1
npy = ny + 1
a = 1.0
b = 0.5
hx = a / nx
hy = b / ny

# Carrega autovetor salvo (φ)
phi = np.loadtxt("phi1.txt")
phi_grid = phi.reshape((npy, npx))  # reshape [linha][coluna]

# Gera malha de coordenadas
x = np.linspace(0, a, npx)
y = np.linspace(0, b, npy)
X, Y = np.meshgrid(x, y)

# ----------- Curvas de nível ------------
plt.figure(figsize=(8, 3))
contour = plt.contourf(X, Y, phi_grid, levels=50, cmap='RdBu_r')
plt.colorbar(contour)
plt.title("Modo φ(x, y)")
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.tight_layout()
plt.savefig("modo_phi1_contour.png")

# ----------- Campo vetorial E = -∇φ ------------
Ex = -np.gradient(phi_grid, axis=1) / hx
Ey = -np.gradient(phi_grid, axis=0) / hy

plt.figure(figsize=(8, 3))
step = 2  # subamostragem para clareza
plt.quiver(X[::step, ::step], Y[::step, ::step], Ex[::step, ::step], Ey[::step, ::step])
plt.title("Campo vetorial E = -∇φ")
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.tight_layout()
plt.savefig("campo_e_vector.png")
