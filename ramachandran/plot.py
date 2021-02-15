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


def get_coordinates_on_reference_map(
        phi_psi_angle: Tuple[float, float],
        reference_map: np.ndarray) -> Tuple[int, int]:

    phi = phi_psi_angle[0]
    psi = phi_psi_angle[1]

    height = reference_map.shape[0]
    width = reference_map.shape[1]

    i = int((180 - psi) / 360 * height)
    j = int((phi + 180) / 360 * width)

    # If i or j == resolution, adjust it.
    if i == height:
        i = height - 1
    if j == width:
        j = width - 1

    return (i, j)


def create_ramachandran_plot(phi_psi_angles: List[Tuple[float, float]],
                             plot_file_path: str,
                             reference_map: Optional[np.ndarray] = None,
                             cmap: Optional[mplcolors.ListedColormap] = None,
                             protein_name: Optional[str] = None,
                             rendering_interpolation: bool = True) -> None:

    phi_psi_angles_numpy = np.array(phi_psi_angles)
    x_numpy = phi_psi_angles_numpy[:, 0]
    y_numpy = phi_psi_angles_numpy[:, 1]

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(111)
    if protein_name is not None:
        ax.set_title(protein_name, fontsize=24)

    interpolation = None
    if rendering_interpolation is True:
        interpolation = "bilinear"

    if reference_map is not None:

        percentile_1 = np.percentile(reference_map, 60)
        percentile_2 = np.percentile(reference_map, 90)

        ax.imshow(np.rot90(reference_map),
                  interpolation=interpolation,
                  cmap=cmap,
                  norm=mplcolors.BoundaryNorm(
                      boundaries=[0, percentile_1, percentile_2, 1],
                      ncolors=cmap.N),
                  origin="upper",
                  extent=(-180, 180, -180, 180))

        # Find outliers
        outliers_idx = []
        for i, phi_psi_angle in enumerate(phi_psi_angles):

            map_i, map_j = get_coordinates_on_reference_map(
                phi_psi_angle=phi_psi_angle,
                reference_map=np.rot90(reference_map))
            if np.rot90(reference_map)[map_i, map_j] < percentile_1:
                outliers_idx.append(i)

        x_outliers_numpy = x_numpy[outliers_idx]
        y_outliers_numpy = y_numpy[outliers_idx]

        x_numpy = np.delete(x_numpy, outliers_idx)
        y_numpy = np.delete(y_numpy, outliers_idx)

        ax.scatter(x_outliers_numpy,
                   y_outliers_numpy,
                   s=20,
                   color="red",
                   edgecolors="black")

    ax.scatter(x_numpy, y_numpy, s=20, color="blue", edgecolors="black")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))

    ax.xaxis.set_tick_params(labelsize=12)
    ax.yaxis.set_tick_params(labelsize=12)

    ax.plot([-180, 180], [0, 0], "--", linewidth=0.5, color="black")
    ax.plot([0, 0], [-180, 180], "--", linewidth=0.5, color="black")

    ax.set_xlabel(r"${\phi}$", fontsize=18, fontweight="bold")
    ax.set_ylabel(r"${\psi}$", fontsize=18, fontweight="bold")

    fig.savefig(plot_file_path, format="svg", dpi=600, bbox_inches="tight")
    plt.close()

    return


def create_ramachandran_plots_from_file(
        file_path: str,
        save_dir_path: str,
        # reference_map_type: Optional[str] = "unsmoothed",
        protein_name: Optional[str] = None,
        rendering_interpolation: bool = False) -> None:

    if not os.path.exists(save_dir_path):
        os.makedirs(save_dir_path)

    residue_torsion_collection = read_residue_torsion_collection_from_file(
        file_path=file_path)

    phi_psi_angles_general = residue_torsion_collection.collect_torsion_angles_general(
    )
    phi_psi_angles_gly = residue_torsion_collection.collect_torsion_angles_gly(
    )
    phi_psi_angles_pro = residue_torsion_collection.collect_torsion_angles_pro(
    )
    phi_psi_angles_prepro = residue_torsion_collection.collect_torsion_angles_prepro(
    )

    phi_psi_angles_list = [
        phi_psi_angles_general, phi_psi_angles_gly, phi_psi_angles_pro,
        phi_psi_angles_prepro
    ]

    package_dir, filename = os.path.split(__file__)

    # Using unsmoothed probability.npz is problematic because
    # many probabilities are exactly zeros and thus the many percentiles are exactly zeros.
    # Plotting these zero values is very problematic.
    # Gaussian density is fine because none of the probability density values are exactly zero.

    # if reference_map_type == "unsmoothed":
    #     npz_file_path = os.path.join(package_dir, "data", "probability.npz")
    #     npz_file = np.load(npz_file_path)
    # elif reference_map_type == "smoothed":
    #     npz_file_path = os.path.join(package_dir, "data", "gaussian_density.npz")
    #     npz_file = np.load(npz_file_path)
    # else:
    #     raise RuntimeError("Unsupported reference map type.")

    npz_file_path = os.path.join(package_dir, "data", "gaussian_density.npz")
    npz_file = np.load(npz_file_path)

    reference_map_general = npz_file["general"]
    reference_map_gly = npz_file["gly"]
    reference_map_pro = npz_file["pro"]
    reference_map_prepro = npz_file["prepro"]

    reference_map_list = [
        reference_map_general, reference_map_gly, reference_map_pro,
        reference_map_prepro
    ]

    # Using Erdős Gábor's cmaps.
    # https://github.com/gerdos/PyRAMA/blob/301df17e5f2c32544b34321c4f8b0254697183ce/pyrama/config.py
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

        create_ramachandran_plot(
            phi_psi_angles=phi_psi_angles,
            reference_map=reference_map,
            cmap=cmap,
            plot_file_path=file_path,
            rendering_interpolation=rendering_interpolation,
            protein_name=protein_name)

    return
