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


def plot_convergence(data_stc, data_eig):
    fig, axs = plt.subplots(3, 2, figsize=(11, 11), sharey="row", sharex="col")
    elems = data_stc["elem"].unique()
    modes = data_eig["mode"].unique()
    data_stc["h"] *= 0.15
    data_eig["h"] *= 0.15

    for elem, ls in zip(elems, ["-", "--"]):
        elem_name = "Triangles" if (elem == 2) else "Quads"
        label = r"{:s}".format(elem_name)
        df = data_stc[data_stc["elem"] == elem]
        kwargs = {"color": f"k", "label": label, "ls": ls, "marker": "o"}
        axs[0, 0].plot(df["h"][1:], df["err_v_ref"][1:], **kwargs)
        axs[0, 1].plot(df["h"], df["err_v_ana"], **kwargs)

    for m in modes:
        for elem, ls in zip(elems, ["-", "--"]):
            elem_name = "Tri" if (elem == 2) else "Quad"
            lb_p = r"$\Phi_{:d}$".format(m) if (elem == 2) else None
            lb_w = r"$\omega_{:d}$".format(m) if (elem == 2) else None
            df = data_eig[(data_eig["elem"] == elem) & (data_eig["mode"] == m)]
            kwargs = {"color": f"C{m-1}", "ls": ls, "marker": "o"}
            axs[1, 0].plot(
                df["h"][1:], df["err_w_ref"][1:], label=lb_w, **kwargs
            )
            axs[2, 0].plot(
                df["h"][1:], df["err_v_ref"][1:], label=lb_p, **kwargs
            )
            axs[1, 1].plot(df["h"], df["err_w_ana"], **kwargs)
            axs[2, 1].plot(df["h"], df["err_v_ana"], **kwargs)

    hs = data_stc["h"].unique()[1:-2]
    label = r"$\mathcal{O}(h^2)$"
    axs[0, 0].plot(hs, hs**2 * 2e-0, ls="-.", label=label, color="gray")
    axs[1, 0].plot(hs, hs**2 * 2e-0, ls="-.", label=label, color="gray")
    axs[2, 0].plot(hs, hs**2 * 5e-1, ls="-.", label=label, color="gray")

    fig.suptitle("FEM convergence - Cantilever beam", fontsize=ftsz_title)
    axs[0, 0].set_title("Error vs fine solution", fontsize=ftsz_axis)
    axs[0, 1].set_title("Error vs 1D beam theory", fontsize=ftsz_axis)

    axs[0, 0].set_ylabel(r"$\|v^h - v^{*}\|$", fontsize=ftsz_axis)
    axs[1, 0].set_ylabel(r"$|\omega^h - \omega^{*}|$", fontsize=ftsz_axis)
    axs[2, 0].set_ylabel(r"$\|\Phi^h - \Phi^{*}\|$", fontsize=ftsz_axis)

    for ax in axs[:, 0]:
        ax.legend(fontsize=ftsz_leg, loc="lower right")

    for ax in axs[-1, :]:
        ax.set_xlabel("$h/H$", fontsize=ftsz_axis)

    for ax in axs.ravel():
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(ls=":", which="both")

    fig.tight_layout()
    # fig.savefig("./tests/conv_beam.pdf", bbox_inches="tight")
    plt.show()
    return


def beam():
    data_stc = load_convergence("./tests/conv_beam_static.txt")
    data_eig = load_convergence("./tests/conv_beam_eigen.txt")
    plot_convergence(data_stc, data_eig)
    return


if __name__ == "__main__":

    beam()
