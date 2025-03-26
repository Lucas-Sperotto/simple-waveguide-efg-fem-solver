import csv
from scipy.special import jn_zeros, jnp_zeros

# Número de ordens m e raízes n
max_m = 10
num_zeros = 10

# Abre o CSV para escrita
with open("bessel_roots.csv", mode="w", newline="") as file:
    writer = csv.writer(file)
    # Cabeçalho
    writer.writerow(["Tipo", "m", "n", "Raiz"])

    for m in range(0, max_m + 1):
        # Raízes para TM (J_m(x) = 0)
        zeros_TM = jn_zeros(m, num_zeros)
        for n, root in enumerate(zeros_TM, start=1):
            writer.writerow(["TM", m, n, f"{root:.15f}"])

        # Raízes para TE (J'_m(x) = 0)
        zeros_TE = jnp_zeros(m, num_zeros)
        for n, root in enumerate(zeros_TE, start=1):
            writer.writerow(["TE", m, n, f"{root:.15f}"])

print("Arquivo 'bessel_roots.csv' gerado com sucesso!")
