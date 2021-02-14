from typing import List, Tuple, Optional
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
import matplotlib.colors as mplcolors
from ramachandran.io import read_residue_torsion_collection_from_file

def create_ramachandran_plot(phi_psi_angles: List[Tuple[float, float]],
                             plot_file_path: str,
                             reference_map: Optional[np.ndarray] = None,
                             cmap: Optional[mplcolors.ListedColormap] = None,
                             protein_name: Optional[str] = None,
                             rendering_interpolation: bool = True) -> None:

    x = []
    y = []

    for phi, psi in phi_psi_angles:
        x.append(phi)
        y.append(psi)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    interpolation = None
    if rendering_interpolation is True:
        interpolation="bilinear"
    
    if reference_map is not None:

        ax.imshow(np.rot90(reference_map),
                interpolation=interpolation,
                cmap=cmap,
                norm=mplcolors.BoundaryNorm(boundaries=[
                    0,
                    np.percentile(reference_map, 50),
                    np.percentile(reference_map, 90), 1
                ],
                                            ncolors=cmap.N),
                origin="upper",
                extent=(-180, 180, -180, 180))

    ax.scatter(x, y, s=10)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))

    ax.plot([-180, 180], [0, 0], "--", linewidth=0.5, color="black")
    ax.plot([0, 0], [-180, 180], "--", linewidth=0.5, color="black")

    ax.set_xlabel(r"${\phi}$", fontsize=16, fontweight="bold")
    ax.set_ylabel(r"${\psi}$", fontsize=16, fontweight="bold")

    fig.savefig(plot_file_path, format="svg", dpi=600, bbox_inches="tight")
    plt.close()

    return


def create_ramachandran_plots_from_file(file_path: str,
                             save_dir_path: str,
                             reference_map_type: Optional[str] = "unsmoothed",
                             protein_name: Optional[str] = None,
                             rendering_interpolation: bool = True) -> None:


    if not os.path.exists(save_dir_path):
        os.makedirs(save_dir_path)

    residue_torsion_collection = read_residue_torsion_collection_from_file(
        file_path=file_path)

    phi_psi_angles_general = residue_torsion_collection.collect_torsion_angles_general()
    phi_psi_angles_gly = residue_torsion_collection.collect_torsion_angles_gly()
    phi_psi_angles_pro = residue_torsion_collection.collect_torsion_angles_pro()
    phi_psi_angles_prepro = residue_torsion_collection.collect_torsion_angles_prepro()

    phi_psi_angles_list = [
        phi_psi_angles_general, phi_psi_angles_gly, phi_psi_angles_pro,
        phi_psi_angles_prepro
    ]

    package_dir, filename = os.path.split(__file__)
    if reference_map_type == "unsmoothed":
        npz_file_path = os.path.join(package_dir, "data", "probability.npz")
        npz_file = np.load(npz_file_path)
    elif reference_map_type == "smoothed":
        npz_file_path = os.path.join(package_dir, "data", "gaussian_density.npz")
        npz_file = np.load(npz_file_path)
    else:
        raise RuntimeError("Unsupported reference map type.")
    
    reference_map_general = npz_file["general"]
    reference_map_gly = npz_file["gly"]
    reference_map_pro = npz_file["pro"]
    reference_map_prepro = npz_file["prepro"]

    reference_map_list = [
        reference_map_general, reference_map_gly, reference_map_pro,
        reference_map_prepro
    ]

    cmap_general = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])
    cmap_gly = mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F'])
    cmap_pro = mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C'])
    cmap_prepro = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])

    cmap_list = [cmap_general, cmap_gly, cmap_pro, cmap_prepro]

    filename_list = ["general.svg", "gly.svg", "pro.svg", "prepro.svg"]
    file_path_list = [
        os.path.join(save_dir_path, filename) for filename in filename_list
    ]

    for phi_psi_angles, reference_map, cmap, file_path in zip(
            phi_psi_angles_list, reference_map_list, cmap_list,
            file_path_list):

        create_ramachandran_plot(phi_psi_angles=phi_psi_angles,
                                 reference_map=reference_map,
                                 cmap=cmap,
                                 plot_file_path=file_path,
                                 rendering_interpolation=rendering_interpolation)


    return 