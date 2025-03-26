import matplotlib.pyplot as plt
import numpy as np
import math

# Parâmetros do guia de onda
a = 1.0
b = 0.5
pi = np.pi

# Gera os primeiros autovalores analíticos TE para m,n ≥ 1
def autovalores_analiticos_TE(a, b, max_index=10):
    valores = []
    for m in range(0, 10):
        for n in range(0, 10):
            if m == 0 and n == 0:
                continue  # Pula o caso (0,0)
            val = math.sqrt((m * pi / a) ** 2 + (n * pi / b) ** 2)
            valores.append((val, f"{m},{n}"))
    valores.sort()
    return valores[:max_index]

aut_analiticos = autovalores_analiticos_TE(a, b, 10)
valores_analiticos = [v[0] for v in aut_analiticos]
rotulos_analiticos = [f"({v[1]})" for v in aut_analiticos]

# Carrega os autovalores FEM e EFGM
fem_vals = np.loadtxt("autovalores_fem.txt")
efgm_vals = np.loadtxt("autovalores_efgm.txt")

# Filtra valores espúrios (≈ 1)
def filtrar(vals, limiar=1.0001):
    return np.array([v for v in vals if v > limiar])

fem_vals_filtrados = np.sort(filtrar(fem_vals))
efgm_vals_filtrados = np.sort(filtrar(efgm_vals))

# Seleciona os primeiros N válidos
N = min(len(fem_vals_filtrados), len(efgm_vals_filtrados), len(valores_analiticos))
x = np.arange(1, N + 1)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x, fem_vals_filtrados[:N], 'o-', label='FEM')
plt.plot(x, efgm_vals_filtrados[:N], 's--', label='EFGM')
plt.plot(x, valores_analiticos[:N], 'd:', label='Analítico TE')

plt.xticks(x, rotulos_analiticos[:N])  # Rótulos (m,n)
plt.xlabel("Modo (m,n)")
plt.ylabel("Autovalores (λ = k²)")
plt.title("Comparação dos Autovalores: FEM x EFGM x Analítico (TE)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("comparacao_autovalores_te_fem_efgm.png")
# plt.show()  # Desnecessário no WSL
