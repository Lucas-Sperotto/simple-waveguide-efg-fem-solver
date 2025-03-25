import matplotlib.pyplot as plt
import numpy as np

# Tamanho da malha (nx elementos → nx+1 nós)
nx = 5
ny = 3
npx = nx + 1
npy = ny + 1

hx = 1.0 / nx
hy = 1.0 / ny

fig, ax = plt.subplots(figsize=(8, 4))

for j in range(npy):
    for i in range(npx):
        x = i * hx
        y = j * hy
        idx = j * npx + i

        # Detecta se é borda
        is_border = (i == 0 or i == npx - 1 or j == 0 or j == npy - 1)
        color = 'red' if is_border else 'black'
        ax.plot(x, y, 'o', color=color)

        # Mostra o índice global do nó
        ax.text(x + 0.01, y + 0.01, f'{idx}', fontsize=9, color=color)

# Grade
for i in range(npx):
    ax.plot([i*hx]*npy, np.linspace(0, 1, npy), '--', color='gray', linewidth=0.5)
for j in range(npy):
    ax.plot(np.linspace(0, 1, npx), [j*hy]*npx, '--', color='gray', linewidth=0.5)

ax.set_aspect('equal')
ax.set_title(f'Nodos da malha estruturada ({nx}x{ny})')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.tight_layout()
plt.grid(True)
plt.savefig("malha_com_bordas.png")
plt.show()
