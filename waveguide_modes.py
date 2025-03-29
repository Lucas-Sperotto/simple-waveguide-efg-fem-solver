import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import meshio
import subprocess
import pandas as pd
import os
from scipy.special import jn_zeros  # Importa função para obter raízes das funções de Bessel

# === Função para gerar a malha circular no GMSH ===
def generate_circle_mesh(radius=0.01, filename='circle'):
    # Cria o conteúdo do arquivo .geo com comandos para GMSH
    geo_content = f"""
    lc = {radius/20};  // Tamanho característico da malha
    Point(1) = {{0, 0, 0, lc}};
    Point(2) = {{{radius}, 0, 0, lc}};
    Point(3) = {{0, {radius}, 0, lc}};
    Point(4) = {{-{radius}, 0, 0, lc}};
    Point(5) = {{0, -{radius}, 0, lc}};
    Circle(1) = {{2, 1, 3}};
    Circle(2) = {{3, 1, 4}};
    Circle(3) = {{4, 1, 5}};
    Circle(4) = {{5, 1, 2}};
    Line Loop(5) = {{1, 2, 3, 4}};
    Plane Surface(6) = {{5}};
    Physical Curve(100) = {{1,2,3,4}};  // Define as curvas da borda
    Physical Surface(200) = {{6}};      // Define a superfície como domínio
    Mesh 2;
    Save "{filename}.msh";              // Salva malha em formato .msh
    """
    geo_file = f"{filename}.geo"
    with open(geo_file, "w") as f:
        f.write(geo_content)  # Escreve o conteúdo no arquivo
    subprocess.run(["gmsh", geo_file, "-2", "-format", "msh2"], check=True)  # Executa o GMSH

# === Função principal para resolver os modos com GMSH e FEM ===
def solve_modes_with_gmsh(radius=0.01, mode='TM', num_modes=6, filename='teste01'):
    # === Geração da malha com GMSH ===
    generate_circle_mesh(radius, filename)
    mesh = meshio.read(f"{filename}.msh")

    # === Leitura dos nós e elementos ===
    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]
    Nn = len(points)  # Número de nós
    X0 = np.zeros(Nn)  # Vetor de solução

    # === Identificação dos nós da borda ===
    boundary_nodes = set()
    for cell_block, phys_ids in zip(mesh.cells, mesh.cell_data_dict["gmsh:physical"].values()):
        if cell_block.type == "line":
            for i, line in enumerate(cell_block.data):
                if phys_ids[i] == 100:
                    boundary_nodes.update(line)

    # === Modo TM: condição de Dirichlet (exclui nós da borda) ===
    if mode == 'TM':
        node_id = np.ones(Nn, dtype=int)
        for i in boundary_nodes:
            node_id[i] = 0
            X0[i] = 0

        # Indexação apenas dos nós internos (não-borda)
        index = np.zeros(Nn, dtype=int)
        counter = 0
        for i in range(Nn):
            if node_id[i] == 1:
                counter += 1
                index[i] = counter
        Nf = counter  # Graus de liberdade

        # Montagem das matrizes S e T com restrição de borda
        S = lil_matrix((Nf, Nf))
        T = lil_matrix((Nf, Nf))
        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)
            b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])
            for i in range(3):
                for j in range(3):
                    Se = (b[i] * b[j] + c_[i] * c_[j]) * Ae
                    Te = Ae / 6 if i == j else Ae / 12
                    if node_id[n[i]] != 0 and node_id[n[j]] != 0:
                        S[index[n[i]]-1, index[n[j]]-1] += Se
                        T[index[n[i]]-1, index[n[j]]-1] += Te

        # Autovalores apenas nos nós internos
        vals, vecs = eigsh(S, k=num_modes+4, M=T, sigma=0, which='LM')
        # === Exporta autovalores, raiz e kc * r ===
        autovalores_completos = np.real(vals)
        raiz_autovalores = np.sqrt(np.clip(autovalores_completos, 0, None))  # Garante não-negatividade
        kc_r = raiz_autovalores * radius

        # Cria DataFrame com autovalores
        df_autovalores = pd.DataFrame({
            'Índice': list(range(1, len(autovalores_completos)+1)),
            'Autovalor (λ)': autovalores_completos,
            'Raiz (√λ = kc)': raiz_autovalores,
            'kc * r': kc_r
        })

        # Salva arquivo em CSV separado
        output_dir = f'out/results/{mode.lower()}_{filename}'
        os.makedirs(output_dir, exist_ok=True)
        autoval_path = os.path.join(output_dir, 'autovalores_detalhados.csv')
        df_autovalores.to_csv(autoval_path, index=False)

        kc = np.sqrt(np.real(vals[:num_modes]))
        fc = 3e8 * kc / (2 * np.pi)

        # Reconstrói a solução completa (incluindo zeros nos nós da borda)
        modos = []
        for q in range(num_modes):
            X0[:] = 0
            j = 0
            for i in range(Nn):
                if index[i] != 0:
                    X0[i] = vecs[j, q]
                    j += 1
            modos.append(X0.copy())

    # === Modo TE: condição de Neumann (todos os nós participam) ===
    elif mode == 'TE':
        # Todos os nós participam da solução (sem exclusão)
        S = lil_matrix((Nn, Nn))
        T = lil_matrix((Nn, Nn))
        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)
            b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])
            for i in range(3):
                for j in range(3):
                    Se = (b[i] * b[j] + c_[i] * c_[j]) * Ae
                    Te = Ae / 6 if i == j else Ae / 12
                    S[n[i], n[j]] += Se
                    T[n[i], n[j]] += Te

        # Autovalores no domínio completo (todos os nós)
        vals, vecs = eigsh(S, k=num_modes+4, M=T, sigma=0, which='LM')
        vals = np.real(vals)

        # === Exporta autovalores, raiz e kc * r ===
        autovalores_completos = np.real(vals)
        raiz_autovalores = np.sqrt(np.clip(autovalores_completos, 0, None))  # Garante não-negatividade
        kc_r = raiz_autovalores * radius

        # Cria DataFrame com autovalores
        df_autovalores = pd.DataFrame({
            'Índice': list(range(1, len(autovalores_completos)+1)),
            'Autovalor (λ)': autovalores_completos,
            'Raiz (√λ = kc)': raiz_autovalores,
            'kc * r': kc_r
        })

        # Salva arquivo em CSV separado
        output_dir = f'out/results/{mode.lower()}_{filename}'
        os.makedirs(output_dir, exist_ok=True)
        autoval_path = os.path.join(output_dir, 'autovalores_detalhados.csv')
        df_autovalores.to_csv(autoval_path, index=False)

        # Elimina o primeiro autovalor (zero, associado ao modo trivial)
        kc = np.sqrt(vals[1:num_modes+1])
        fc = 3e8 * kc / (2 * np.pi)

        # Coleta os modos correspondentes
        modos = [vecs[:, i] for i in range(1, num_modes+1)]

    else:
        raise ValueError("Modo inválido. Use 'TM' ou 'TE'.")

    # === Diretório para salvar as imagens ===
    save_path = f"out/img/{mode.lower()}_{filename}"
    os.makedirs(save_path, exist_ok=True)

    # === Geração dos gráficos ===
    for q in range(num_modes):
        plt.figure()
        plt.tricontourf(points[:, 0], points[:, 1], triangles, modos[q], levels=100, cmap='jet')
        plt.colorbar()
        plt.axis('equal')
        plt.tight_layout()
        plt.title(f"Modo {q+1} - {mode} | fc = {fc[q]/1e9:.3f} GHz")
        fig_path = os.path.join(save_path, f"modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()

    return fc, modos


# === Exporta os dados em CSV ===
def export_results_to_csv(fc, fcreal, error, kc, radius, mode='TM', filename='teste01'):
    # Calcula kc * r (valor adimensional)
    kc_r = kc * radius
    
    # Monta o dicionário com os dados
    data = {
        'Modo': [f'{i+1}' for i in range(len(fc))],
        'Frequência FEM (GHz)': fc / 1e9,
        'Frequência Teórica (GHz)': fcreal / 1e9,
        'kc * r': kc_r,
        'Erro Relativo (%)': error
    }

    # Cria DataFrame e salva em CSV
    df = pd.DataFrame(data)
    output_dir = f'out/results/{mode.lower()}_{filename}'
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'frequencias_modos.csv')
    df.to_csv(csv_path, index=False)
    return df


# === Bloco principal de execução ===
from scipy.special import jn_zeros

if __name__ == "__main__":
    # Parâmetros
    num_modos = 24
    raio = 0.01  # 10 mm

    # === Cálculo para modos TE ===
    print("Calculando modos TE...")
    fc_te, modos_te = solve_modes_with_gmsh(radius=raio, mode='TE', num_modes=num_modos, filename='te_24modos')

    # Teóricos para TE: raízes de J1
    pnm_TE = jn_zeros(1, num_modos)
    fcreal_TE = 3e8 * pnm_TE / (2 * np.pi * raio)
    kc_te = 2 * np.pi * fcreal_TE / 3e8  # kc = 2πf / c
    kc_te_fem = 2 * np.pi * fc_te / 3e8
    erro_TE = 100 * np.abs((fcreal_TE - fc_te)) / fcreal_TE

    # Exportar resultados TE
    df_te = export_results_to_csv(fc_te, fcreal_TE, erro_TE, kc_te_fem, raio, mode='TE', filename='te_24modos')
    print("Resultados dos modos TE:")
    print(df_te)

    # === Cálculo para modos TM ===
    print("\nCalculando modos TM...")
    fc_tm, modos_tm = solve_modes_with_gmsh(radius=raio, mode='TM', num_modes=num_modos, filename='tm_24modos')

    # Teóricos para TM: raízes de J0
    pnm_TM = jn_zeros(0, num_modos)
    fcreal_TM = 3e8 * pnm_TM / (2 * np.pi * raio)
    kc_tm = 2 * np.pi * fcreal_TM / 3e8
    kc_tm_fem = 2 * np.pi * fc_tm / 3e8
    erro_TM = 100 * np.abs((fcreal_TM - fc_tm)) / fcreal_TM

    # Exportar resultados TM
    df_tm = export_results_to_csv(fc_tm, fcreal_TM, erro_TM, kc_tm_fem, raio, mode='TM', filename='tm_24modos')
    print("Resultados dos modos TM:")
    print(df_tm)
