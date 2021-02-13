from typing import List, Tuple, Optional
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
import matplotlib.colors as mplcolors
from ramachandran.read_structure import read_pdb, read_pdbx
from ramachandran.compute_dihedral_angle import protein_backbone_dihedral_angle_phi, protein_backbone_dihedral_angle_psi, collect_dihedral_angles

def create_ramachandran_trusted_regsion_plot(trusted_regsion_filepath: str, plot_filepath: str) -> None:

    npzfile = np.load(trusted_regsion_filepath)

    counts_density = npzfile["density_map"]

    print(np.sum(counts_density))


    cmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.imshow(np.rot90(counts_density), cmap=cmap, norm=mplcolors.BoundaryNorm([0, 5e-7, 2e-5, 1], cmap.N), origin="upper", extent=(-180, 180, -180, 180))

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))

    ax.plot([-180, 180], [0, 0], "--", linewidth=0.5, color="black")
    ax.plot([0, 0], [-180, 180], "--", linewidth=0.5, color="black")

    ax.set_xlabel(r"${\phi}$", fontsize=16, fontweight="bold")
    ax.set_ylabel(r"${\psi}$", fontsize=16, fontweight="bold")

    fig.savefig(plot_filepath, format="svg", dpi=600, bbox_inches="tight")
    plt.close()


    return



def create_ramachandran_plot(filepath: str,
                             plot_filepath: str,
                             protein_name: Optional[str] = None) -> None:

    dihedral_angles = collect_dihedral_angles(filepath=filepath)

    x = []
    y = []

    for phi, psi in dihedral_angles:
        x.append(phi)
        y.append(psi)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=10)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))

    ax.plot([-180, 180], [0, 0], "--", linewidth=0.5, color="black")
    ax.plot([0, 0], [-180, 180], "--", linewidth=0.5, color="black")

    ax.set_xlabel(r"${\phi}$", fontsize=16, fontweight="bold")
    ax.set_ylabel(r"${\psi}$", fontsize=16, fontweight="bold")

    fig.savefig(plot_filepath, format="svg", dpi=600, bbox_inches="tight")
    plt.close()

    return
