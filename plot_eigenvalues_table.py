import matplotlib.pyplot as plt
import numpy as np

# Parâmetros do guia de onda
a = 1.0
b = 0.5
pi = np.pi

# Gera os primeiros autovalores analíticos TE para m,n ≥ 1
def autovalores_analiticos_TE(a, b, max_index=10):
    valores = []
    for m in range(1, 10):
        for n in range(1, 10):
            val = (m * pi / a) ** 2 + (n * pi / b) ** 2
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
def filtrar(vals, limiar=1.01):
    return np.array([v for v in vals if v > limiar])

fem_vals_filtrados = np.sort(filtrar(fem_vals))
efgm_vals_filtrados = np.sort(filtrar(efgm_vals))

# Seleciona os primeiros N válidos
N = min(len(fem_vals_filtrados), len(efgm_vals_filtrados), len(valores_analiticos))

# Calcular erro relativo (%)
def erro_relativo(calc, ref):
    return np.abs((calc - ref) / ref) * 100

# Tabela com os dados
tabela_dados = []
for i in range(N):
    analitico = valores_analiticos[i]
    fem = fem_vals_filtrados[i]
    efgm = efgm_vals_filtrados[i]
    erro_fem = erro_relativo(fem, analitico)
    erro_efgm = erro_relativo(efgm, analitico)
    tabela_dados.append([
        rotulos_analiticos[i],
        f"{analitico:.4f}",
        f"{fem:.4f}",
        f"{erro_fem:.2f}%",
        f"{efgm:.4f}",
        f"{erro_efgm:.2f}%"
    ])

# Cabeçalho da tabela
col_labels = ["Modo (m,n)", "Analítico", "FEM", "Erro FEM", "EFGM", "Erro EFGM"]

# Criar figura da tabela
fig, ax = plt.subplots(figsize=(12, 0.5 * N + 1))
ax.axis('off')
tabela = ax.table(cellText=tabela_dados, colLabels=col_labels, loc='center', cellLoc='center')
tabela.auto_set_font_size(False)
tabela.set_fontsize(10)
tabela.scale(1.2, 1.2)

plt.title("Tabela Comparativa dos Autovalores e Erros Relativos (TE)", pad=20)
plt.tight_layout()
plt.savefig("tabela_autovalores_te.png")
# plt.show()  # Desnecessário em ambiente WSL
