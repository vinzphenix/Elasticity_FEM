"""
File : perf_CG.py
Author : Vincent Degrooff
Date : 2025

Description :
    Display results of the CG complexity analysis
    
Project :
    FEM Simulation Toolkit for Linear Elasticity
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ftsz_axis = 16
ftsz_title = 18
plt.rcParams["font.family"] = "monospace"


def load_data(filename):
    data = pd.read_csv(filename + ".txt", delim_whitespace=True)
    data["time"] = data["tprec"] + data["tsolve"]
    return data


def plot_data(data, filename):
    dim = int(filename[-2])
    n_cut = 1e3 if dim == 1 else 1e5
    data = data[data["n"] >= n_cut]
    x_cplx = np.linspace(data["n"].min(), data["n"].max(), 5)
    x_ticks = np.array([1e5, 2e5, 5e5, 1e6])
    y_ticks = np.array([1.0, 2.0, 5.0, 10.0])
    fmt = ".0f"
    
    if dim == 1:
        suptitle = "Regular periodic line grid"
        complexity = 1e-9 * x_cplx[1:] ** 2
        label = "$O(n^2)$"
        n_to_nx = lambda n: n
        bounds = (0.9e3, 4.e4)
        x_ticks = 2.5 * x_ticks // 1e2
        y_ticks = np.power(10., np.arange(-4, 2))
        fmt = "g"
    elif dim == 2:
        suptitle = "Regular periodic quad grid"
        complexity = 1e-8 * x_cplx[1:] ** 1.5
        label = "$O(n^{1.5})$"
        n_to_nx = lambda n: np.sqrt(n).astype(int)
        bounds = (0.95e5, 1.1e6)
    elif dim == 3:
        suptitle = "Regular periodic hex grid"
        complexity = 2e-7 * x_cplx[1:] ** 1.25
        label = "$O(n^{1.25})$"
        n_to_nx = lambda n: np.cbrt(n).astype(int)
        bounds = (0.95e5, 1.1e6)
    else:
        raise ValueError("Invalid dimension")

    fig, axs = plt.subplots(1, 3, figsize=(16, 6))
    methods = data["method"].unique()
    for method in methods:
        lw = 2.5 if method == "NoPrec" else 1.5
        df = data[data["method"] == method]
        axs[0].plot(df["n"], df["time"], "-o", lw=lw, label=f"{method}")
        axs[1].plot(df["n"], df["it"], "-o", lw=lw, label=f"{method}")
        ratio_prec = 1e2 * df["tprec"] / df["time"]
        axs[2].plot(df["n"], ratio_prec, "-o", lw=lw, label=f"{method}")

    axs[0].plot(x_cplx[1:], complexity, ls="--", label=label, color="gray")
    axs[0].set_yscale("log")
    axs[0].set_yticks(y_ticks)
    axs[0].set_yticklabels([f"{yt:{fmt}}" for yt in y_ticks])

    axs[0].set_title("CG total runtime [s]", fontsize=ftsz_axis)
    axs[1].set_title("CG iterations", fontsize=ftsz_axis)
    axs[2].set_title("Preconditioner runtime [%]", fontsize=ftsz_axis)
    fig.suptitle(suptitle, fontsize=ftsz_title)

    x_tick_labels = [f"{xt:.0e}" for xt in x_ticks]
    h_ticks = 1.0 / n_to_nx(x_ticks)
    for ax in axs:
        ax.set_xscale("log")
        ax.set_xlim(*bounds)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels)
        ax.set_xlabel("n", fontsize=ftsz_axis)
        color = "grey"
        ax2 = ax.twiny()
        ax2.set_xscale("log")
        ax2.set_xbound(ax.get_xlim())
        ax2.set_xticks(x_ticks, minor=False, color="C1")
        ax2.tick_params(axis="x", colors=color)
        ax2.set_xticks([], minor=True)
        ax2.set_xticklabels([f"{h:.1e}" for h in h_ticks], color=color)
        # ax2.set_xlabel(r"Step size h", fontsize=ftsz_axis, color=color)
        ax.grid(ls=":", which="both")
        ax.legend()

    fig.tight_layout()
    # fig.savefig(filename + ".svg", bbox_inches="tight")
    plt.show()
    return


if __name__ == "__main__":
    filename = "./perf_CG_3D"
    data = load_data(filename)
    plot_data(data, filename)
