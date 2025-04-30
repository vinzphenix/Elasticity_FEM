"""
File : conv_tank.py
Author : Vincent Degrooff
Date : 2025

Description :
    Display results of the thick-wall cylinder convergence analysis
    
Project :
    FEM Simulation Toolkit for Linear Elasticity
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ftsz_leg = 14
ftsz_axis = 16
ftsz_title = 20
# plt.rcParams["font.family"] = "monospace"
plt.rcParams["text.usetex"] = True


def load_convergence(filename):
    data = pd.read_csv(filename, delim_whitespace=True)
    return data


def plot_convergence(data):
    fig, axs = plt.subplots(3, 2, figsize=(11, 11), sharey="row", sharex="col")
    elems = data["etype"].unique()
    modes = data["mtype"].unique()
    data["h"] *= 0.014 / 0.05

    for i, mode in enumerate(modes):
        for elem, ls in zip(elems, ["-", "--"]):
            elem_name = "Triangles" if (elem == 2) else "Quads"
            label = r"{:s}".format(elem_name)
            df = data[(data["mtype"] == mode) & (data["etype"] == elem)]
            kwargs = {"color": f"k", "label": label, "ls": ls, "marker": "o"}
            axs[0, i].plot(df["h"], df["u_r"], **kwargs)
            kwargs = {"ls": ls, "marker": "o"}
            lb_rr = (r"$\sigma_{rr}$" if (elem == 2) else None)
            lb_tt = (r'$\sigma_{\theta\theta}$' if (elem == 2) else None)
            lb_zz = (r"$\sigma_{zz}$" if (elem == 2) else None)
            axs[1, i].plot(df["h"], df["s_lsq_rr"], 'C0', label=lb_rr, **kwargs)
            axs[1, i].plot(df["h"], df["s_lsq_tt"], 'C1', label=lb_tt, **kwargs)
            axs[1, i].plot(df["h"], df["s_lsq_zz"], 'C2', label=lb_zz, **kwargs)
            axs[2, i].plot(df["h"], df["s_avg_rr"], 'C0', label=lb_rr, **kwargs)
            axs[2, i].plot(df["h"], df["s_avg_tt"], 'C1', label=lb_tt, **kwargs)
            axs[2, i].plot(df["h"], df["s_avg_zz"], 'C2', label=lb_zz, **kwargs)

    hs = data["h"].unique()[:-2]
    for ax in axs[0, :]:
        label = r"$\mathcal{O}(h^2)$"
        ax.plot(hs, hs**2 * 2e-2, ls="-.", label=label, color="gray")
    for ax in axs[1:, :].ravel():
        label = r"$\mathcal{O}(h^{4/3})$"
        ax.plot(hs, hs**(4./3.) * 5e-1, ls="-.", label=label, color="gray")
        
    fig.suptitle("FEM convergence - Thick wall cylinder", fontsize=ftsz_title)
    kwargs = {"fontsize": ftsz_axis}
    axs[0, 0].set_title("2D Plane strain", **kwargs)
    axs[0, 1].set_title("3D Axisymmetric", **kwargs)
    axs[0, 0].set_ylabel(r"$\|   u^h -    u^{*}\|_{L_2}$", **kwargs)
    axs[1, 0].set_ylabel(r"$\|\sigma^{\textrm{lsq}} - \sigma^{*}\|_{L_2}$", **kwargs)
    axs[2, 0].set_ylabel(r"$\|\sigma^{\textrm{avg}} - \sigma^{*}\|_{L_2}$", **kwargs)

    for ax in axs[:, 0]:
        ax.legend(fontsize=ftsz_leg, loc="lower right")
        
    for ax in axs[-1, :]:
        ax.set_xlabel(r"$h/\Delta r$", fontsize=ftsz_axis)

    for ax in axs.ravel():
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(ls=":", which="both")

    fig.tight_layout()
    # fig.savefig("./analysis/conv_tank.pdf", bbox_inches="tight")
    plt.show()
    return


def tank():
    data = load_convergence("./analysis/conv_tank.txt")
    plot_convergence(data)
    return


if __name__ == "__main__":

    tank()
