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
        plt.figure(figsize=(6, 5))
        plt.tricontourf(points[:, 0], points[:, 1], triangles, modos[q], levels=100, cmap='jet')
        plt.colorbar(label='Amplitude')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.axis('equal')
        plt.title(f"Modo {q+1} - {mode}\nfc = {fc[q]/1e9:.3f} GHz | kc·r = {kc[q]*radius:.3f}")
#plt.title(f"Modo {q+1} - {mode}\nfc = {fc[q]/1e9:.3f} GHz | kc·r = {kc[q] * radius:.3f}")
        plt.tight_layout()
        fig_path = os.path.join(save_path, f"modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()





    # === Pasta para gráficos vetoriais ===
    save_path_quiver = f"out/img/quiver_{mode.lower()}_{filename}"
    os.makedirs(save_path_quiver, exist_ok=True)

    # === Calcula gradiente e gera quiver plot para cada modo ===
    from scipy.spatial import cKDTree

    tree = cKDTree(points)
    delta = radius / 1000  # Passo pequeno para estimar derivadas

    # === Novo cálculo vetorial usando interpolação FEM (gradiente por triângulo) ===
    save_path_quiver = f"out/img/quiver_{mode.lower()}_{filename}"
    os.makedirs(save_path_quiver, exist_ok=True)

    from scipy.constants import mu_0, pi

    for q in range(num_modes):
        campo_z = modos[q]
        Ex = np.zeros(Nn)
        Ey = np.zeros(Nn)
        contagem = np.zeros(Nn)  # Contar quantas vezes cada nó participa (para média)

        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)

            b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])

            grad_phi_x = np.dot(campo_z[n], b)  # ∂φ/∂x
            grad_phi_y = np.dot(campo_z[n], c_)  # ∂φ/∂y

            if mode == 'TM':
                Ex_local = -grad_phi_y
                Ey_local = grad_phi_x
            elif mode == 'TE':
                kc_local = kc[q] if kc[q] != 0 else 1e-12
                omega = 2 * pi * fc[q]
                Z = omega * mu_0 / kc_local
                Ex_local = -Z * grad_phi_y
                Ey_local = -Z * grad_phi_x

            for i in n:
                Ex[i] += Ex_local
                Ey[i] += Ey_local
                contagem[i] += 1

        # Média dos valores nos nós
        Ex /= np.maximum(contagem, 1)
        Ey /= np.maximum(contagem, 1)

       
                # Cálculo da magnitude para normalização
        magnitude = np.sqrt(Ex**2 + Ey**2)
        nonzero = magnitude > 1e-14
        Ex_unit = np.zeros_like(Ex)
        Ey_unit = np.zeros_like(Ey)
        Ex_unit[nonzero] = Ex[nonzero] / magnitude[nonzero]
        Ey_unit[nonzero] = Ey[nonzero] / magnitude[nonzero]

        # Coloração pela magnitude original
        color_data = magnitude
        

                # === Cálculo da magnitude original dos vetores ===
        magnitude = np.sqrt(Ex**2 + Ey**2)
        max_mag = np.max(magnitude)

        # === Normalização de direção dos vetores ===
        Ex_unit = np.zeros_like(Ex)
        Ey_unit = np.zeros_like(Ey)
        nonzero = magnitude > 1e-14
        Ex_unit[nonzero] = Ex[nonzero] / magnitude[nonzero]
        Ey_unit[nonzero] = Ey[nonzero] / magnitude[nonzero]

        # === Ajuste de magnitude para visualização ===
        magnitude_norm = np.zeros_like(magnitude)
        magnitude_norm[nonzero] = magnitude[nonzero] / max_mag
        magnitude_scaled = 0.3 + 0.7 * magnitude_norm  # mínimo 30% do vetor máximo

        # === Vetores finais com tamanho ajustado ===
        Ex_final = Ex_unit * magnitude_scaled
        Ey_final = Ey_unit * magnitude_scaled

        # === Quiver Plot ===
        plt.figure(figsize=(6, 5))
        quiv = plt.quiver(
            points[:, 0], points[:, 1],
            Ex_final, Ey_final,
            magnitude,  # color map by original magnitude
            cmap='viridis',
            scale=30, width=0.002, edgecolor='k', linewidth=0.1
        )
        plt.colorbar(quiv, label='|Eₜ| (não normalizado)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        # Contorno do guia
        theta = np.linspace(0, 2 * np.pi, 200)
        x_circ = radius * np.cos(theta)
        y_circ = radius * np.sin(theta)
        plt.plot(x_circ, y_circ, 'k--', linewidth=1)

        plt.title(f"Campo Transversal - Modo {q+1} - {mode}\n"
                  f"fc = {fc[q]/1e9:.3f} GHz | kc·r = {kc[q]*radius:.3f}")
        plt.axis('equal')
        plt.tight_layout()
        fig_path = os.path.join(save_path_quiver, f"quiver_modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
     
       
       
       
       
       
       
       
        # === Quiver plot ===
        #plt.figure(figsize=(6, 5))
        # Escalonamento automático baseado na magnitude
        #magnitude = np.sqrt(Ex**2 + Ey**2)
        #max_magnitude = np.max(magnitude)
        #if max_magnitude > 0:
        #    Ex_scaled = Ex / max_magnitude
        #    Ey_scaled = Ey / max_magnitude
        #else:
        #    Ex_scaled = Ex
       #     Ey_scaled = Ey

        #scale_factor = 20  # Você pode ajustar esse valor conforme o visual desejado
        #plt.quiver(points[:, 0], points[:, 1], Ex_scaled, Ey_scaled, scale=scale_factor, width=0.002)
        #plt.xlabel('x (m)')
        #plt.ylabel('y (m)')
        #plt.title(f"Campo Transversal - Modo {q+1} - {mode}")
        #plt.axis('equal')
        #plt.tight_layout()
        #fig_path = os.path.join(save_path_quiver, f"quiver_modo_{q+1}_{mode}.png")
        #plt.savefig(fig_path, dpi=300)
        #plt.close()







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
