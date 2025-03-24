import matplotlib.pyplot as plt
import numpy as np

# Carrega os autovalores
fem_vals = np.loadtxt("autovalores_fem.txt")
efgm_vals = np.loadtxt("autovalores_efgm.txt")

# Ordena os valores
fem_vals.sort()
efgm_vals.sort()

# Seleciona os primeiros N autovalores
N = min(len(fem_vals), len(efgm_vals), 10)
x = np.arange(1, N+1)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(x, fem_vals[:N], 'o-', label='FEM')
plt.plot(x, efgm_vals[:N], 's--', label='EFGM')
plt.xlabel("Índice do modo")
plt.ylabel("Autovalores (λ = k²)")
plt.title("Comparação dos Autovalores FEM x EFGM")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("comparacao_autovalores.png")
plt.show()
