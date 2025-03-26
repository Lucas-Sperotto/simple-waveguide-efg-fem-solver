import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import meshio
import subprocess
import pandas as pd
import os

def generate_circle_mesh(radius=0.01, filename='circle'):
    geo_content = f"""
    lc = {radius/20};
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
    Physical Curve(100) = {{1,2,3,4}};
    Physical Surface(200) = {{6}};
    Mesh 2;
    Save "{filename}.msh";
    """
    geo_file = f"{filename}.geo"
    with open(geo_file, "w") as f:
        f.write(geo_content)
    subprocess.run(["gmsh", geo_file, "-2", "-format", "msh2"], check=True)

def solve_modes_with_gmsh(radius=0.01, mode='TM', num_modes=6, filename='teste01'):
    generate_circle_mesh(radius, filename)
    mesh = meshio.read(f"{filename}.msh")

    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]
    Nn = len(points)
    Ne = len(triangles)
    X0 = np.zeros(Nn)
    node_id = np.ones(Nn, dtype=int)

    # === Encontrar nós da borda ===
    boundary_nodes = set()
    for cell_block, phys_ids in zip(mesh.cells, mesh.cell_data_dict["gmsh:physical"].values()):
        if cell_block.type == "line":
            for i, line in enumerate(cell_block.data):
                if phys_ids[i] == 100:
                    boundary_nodes.update(line)
    for i in boundary_nodes:
        node_id[i] = 0
        X0[i] = 0

    # === Indexação dos nós desconhecidos ===
    index = np.zeros(Nn, dtype=int)
    counter = 0
    for i in range(Nn):
        if node_id[i] == 1:
            counter += 1
            index[i] = counter
    Nf = counter

    # === Montagem das matrizes ===
    S = lil_matrix((Nf, Nf))
    T = lil_matrix((Nf, Nf))
    for tri in triangles:
        n = tri
        x1, y1 = points[n[0]]
        x2, y2 = points[n[1]]
        x3, y3 = points[n[2]]
        mat = np.array([[1, x1, y1], [1, x2, y2], [1, x3, y3]])
        De = np.linalg.det(mat)
        Ae = abs(De / 2)
        b = np.array([(y2 - y3) / De, (y3 - y1) / De, (y1 - y2) / De])
        c_ = np.array([(x3 - x2) / De, (x1 - x3) / De, (x2 - x1) / De])
        for i in range(3):
            for j in range(3):
                Se = (b[i] * b[j] + c_[i] * c_[j]) * Ae
                Te = Ae / 6 if i == j else Ae / 12
                if node_id[n[i]] != 0 and node_id[n[j]] != 0:
                    S[index[n[i]]-1, index[n[j]]-1] += Se
                    T[index[n[i]]-1, index[n[j]]-1] += Te

    # === Resolução dos autovalores ===
    vals, vecs = eigsh(S, k=num_modes+4, M=T, sigma=0, which='LM')
    kc = np.sqrt(np.real(vals[:num_modes]))
    fc = 3e8 * kc / (2 * np.pi)

    # === Teórico ===
    if mode == 'TM':
        pnm = np.array([7.016, 6.38, 5.52, 5.135, 3.832, 2.405])
    else:
        pnm = np.array([5.317, 4.201, 3.832, 3.054, 1.841, 0])
    fcreal = 3e8 * pnm / (2 * np.pi * radius)
    error = 100 * np.abs((fcreal - fc)) / fcreal

    # === Diretório de saída ===
    save_path = f"out/img/{mode.lower()}_{filename}"
    os.makedirs(save_path, exist_ok=True)

    # === Geração das figuras ===
    for q in range(num_modes):
        X0[:] = 0
        j = 0
        for i in range(Nn):
            if index[i] != 0:
                X0[i] = vecs[j, q]
                j += 1

        plt.figure()
        plt.tricontourf(points[:, 0], points[:, 1], triangles, X0, levels=100, cmap='jet')
        plt.colorbar()
        plt.axis('equal')
        plt.tight_layout()

        # Título com todas as infos
        fnum = fc[q] / 1e9
        ftheo = fcreal[q] / 1e9
        err = error[q]
        plt.title(f"Modo {q+1} - {mode} | fc = {fnum:.3f} GHz | fₜₕ = {ftheo:.3f} GHz | erro = {err:.2f} %")

        # Caminho do arquivo
        fig_path = os.path.join(save_path, f"modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()

    return fc, error

def export_results_to_csv(fc, fcreal, error, mode='TM', filename='teste01'):
    data = {
        'Modo': [f'{i+1}' for i in range(len(fc))],
        'Frequência FEM (GHz)': fc / 1e9,
        'Frequência Teórica (GHz)': fcreal / 1e9,
        'Erro Relativo (%)': error
    }
    df = pd.DataFrame(data)
    output_dir = f'out/results/{mode.lower()}_{filename}'
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'frequencias_modos.csv')
    df.to_csv(csv_path, index=False)
    return df

# === Execução adicional para exportar os resultados ===
if __name__ == "__main__":
    fc, error = solve_modes_with_gmsh(radius=0.01, mode='TM', num_modes=6, filename='teste01')
    pnm_TM = np.array([7.016, 6.38, 5.52, 5.135, 3.832, 2.405])
    fcreal = 3e8 * pnm_TM / (2 * np.pi * 0.01)
    error = 100 * np.abs((fcreal - fc)) / fcreal
    df = export_results_to_csv(fc, fcreal, error, mode='TM', filename='teste01')
    print(df)