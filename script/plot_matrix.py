import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable

np.set_printoptions(
    precision=4, linewidth=200, formatter={"float": "{: 8.4e}".format}, edgeitems=10
)


def plot_matrix(systems, titles):
    fig, axs = plt.subplots(
        1, len(systems), figsize=(10, 5), sharex=True, sharey=True
    )
    fig.canvas.manager.set_window_title("Matrices sparsity pattern")
    if type(axs) is not np.ndarray:
        axs = [axs]
    vmax = 0.
    for matrix, _, mmax in systems:
        vmax = max(vmax, mmax)
    for ax, (matrix, _, mmax), title in zip(axs, systems, titles):
        ax.set_title(title)
        cs = ax.imshow(matrix.toarray(), cmap="bwr", vmin=-vmax, vmax=vmax)
        # ax.spy(matrix, markersize=2)
        ax.set_aspect("equal")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", pad=0.05, size=0.25)
    fig.colorbar(cs, cax=cax)
    fig.tight_layout()
    # ax.axis("off")
    plt.show()


def compute_eigenvalues(A, M, title):
    print(title)
    n_eigws = 9
    if sp.issparse(A):
        eigws, eigvs = sp.linalg.eigsh(A, k=n_eigws, M=M, which="LM")
        print("Large eigenvalues : ", eigws)
        eigws, eigvs = sp.linalg.eigsh(A, k=n_eigws, M=M, which="SM")
        print("Small eigenvalues : ", eigws)
    else:
        eigws, eigvs = np.linalg.eigh(A)
        print("Eigenvalues : ", eigws)
    print()
    return


def load_matrix(file_name):
    rhs = None
    with open(file_name, "r") as f:
        # Read the first line : # kind n k
        first = f.readline().strip("# ")
        args = first.split()
        kind = args[0]
        if kind == "FULL":
            n = int(args[1])
            rhs_flag = int(args[2])
            matrix = np.loadtxt(f)
            if rhs_flag == 1:
                rhs = matrix[:, -1]
                matrix = matrix[:, :-1]
            max_abs = np.amax(np.abs(matrix))
        elif kind == "SYM_BAND":
            n = int(args[1])
            k = int(args[2])
            rhs_flag = int(args[3])
            data = np.loadtxt(f)
            diags = [data[k - i : n, i] for i in range(k + 1)]
            diags += [data[i:n, k - i] for i in range(1, k + 1)]
            matrix = sp.diags(
                diags, offsets=np.arange(-k, k + 1), shape=(n, n)
            )
            if rhs_flag == 1:
                rhs = data[:, -1]
            max_abs = np.amax(np.abs(matrix).data)
        elif kind == "CSR":
            n = int(args[1])
            nnz = int(args[2])
            rhs_flag = int(args[3])
            memory = np.loadtxt(f, max_rows=n).reshape((n, -1))
            row_ptr = np.empty(n + 1, dtype=int)
            row_ptr[:-1] = memory[:, 0].astype(int)
            row_ptr[-1] = nnz
            if rhs_flag == 1:
                rhs = memory[:, -1]
            col_data = np.loadtxt(f, max_rows=nnz)
            cols = col_data[:, 0].astype(int)
            data = col_data[:, 1]
            matrix = sp.csr_matrix((data, cols, row_ptr), shape=(n, n))
            matrix = matrix + matrix.T - sp.diags(matrix.diagonal())
            sum_lines = np.array(np.abs(matrix).sum(axis=1)).ravel()
            print("Sum_j |K_ij|   |   Diagonal")
            print(np.c_[sum_lines, matrix.diagonal()])
            max_abs = np.amax(np.abs(matrix).data)
            
    return matrix, rhs, max_abs


if __name__ == "__main__":
    # mass = load_matrix("../M.mtx")
    # compute_eigenvalues(mass, None, "Mass Matrix")
    # stiff = load_matrix("../K.mtx")
    # compute_eigenvalues(stiff, None, "Stiffness Matrix")
    # compute_eigenvalues(stiff, mass, "Generalized Eigenvalues")
    # matrices = [mass, stiff]
    # names = ["Mass Matrix", "Stiffness Matrix"]

    # stiff = load_matrix("./K_csr.txt")
    stiff = load_matrix("./K.txt")
    compute_eigenvalues(stiff[0].toarray(), None, "Stiffness Matrix")
    
    matrix = stiff[0].toarray()
    sol = np.linalg.solve(matrix, stiff[1])
    print("Solution : ")
    print(sol)
    
    
    matrices = [stiff]
    names = ["Stiffness Matrix"]
    plot_matrix(matrices, names)
