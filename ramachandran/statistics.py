from typing import List, Tuple, Optional
import os
import math
from tqdm import tqdm
from multiprocessing import Pool
import numpy as np
from scipy import stats
from ramachandran.io import read_residue_torsion_collection_from_directory


def count_torsion_angles(phi_psi_angles: List[Tuple[float, float]],
                         resolution: int = 90,
                         normalized: bool = False) -> np.ndarray:
    """Count the torsion angles for histogram. No smoothing function was applied.

    Args:
        phi_psi_angles (List[Tuple[float, float]]): Torsion phi psi angles.
        resolution (int, optional): Resolution of the counts. Defaults to 90.
        normalized (bool, optional): Normalize the counts to density. Defaults to False.

    Returns:
        np.ndarray: A Numpy array of shape [resolution, resolution].
    """

    counts = np.zeros(shape=(resolution, resolution), dtype=np.uint32)

    for phi, psi in phi_psi_angles:

        i = int((phi + 180) / 360 * resolution)
        j = int((psi + 180) / 360 * resolution)
        # i = int((180 - psi) / 360 * resolution)
        # j = int((phi + 180) / 360 * resolution)

        # If i or j == resolution, adjust it.
        if i == resolution:
            i = resolution - 1
        if j == resolution:
            j = resolution - 1

        counts[i, j] += 1

    if normalized == True:
        counts = counts.astype(np.float32) / np.sum(counts)

    return counts


def compute_gaussian_kde_density(phi_psi_angles: List[Tuple[float, float]],
                                 resolution: int = 180) -> np.ndarray:
    """Use Gaussian kernel density smoothed estimation for the torsion angles.

    Args:
        phi_psi_angles (List[Tuple[float, float]]): Torsion phi psi angles.
        resolution (int, optional): Resolution of the density. Defaults to 180.

    Returns:
        np.ndarray: A Numpy array of shape [resolution, resolution].
    """

    # Turn the list of tuples to Numpy array
    # Shape: [2, len(phi_psi_angles)]
    values = np.array(phi_psi_angles).T

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    xmin = -180
    xmax = 180
    ymin = -180
    ymax = 180

    X, Y = np.mgrid[xmin:xmax:(1j * resolution), ymin:ymax:(1j * resolution)]
    positions = np.vstack([X.ravel(), Y.ravel()])

    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    return Z


def count_torsion_angles_from_directory(
    dir_path: str,
    save_file_path: Optional[str] = None,
    b_factor_threshold: float = 30,
    resolution: int = 180,
    num_processes: int = 4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    residue_torsion_collection = read_residue_torsion_collection_from_directory(
        dir_path=dir_path,
        b_factor_threshold=b_factor_threshold,
        num_processes=num_processes)

    # Counting residue torsions
    # GLY
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_gly()
    print("Couting torsions for GLY, {} samples, this might take a while...".
          format(len(phi_psi_angles)))
    counts_gly = count_torsion_angles(phi_psi_angles=phi_psi_angles,
                                      resolution=resolution,
                                      normalized=True)
    # PRO
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_pro()
    print("Couting torsions for PRO, {} samples, this might take a while...".
          format(len(phi_psi_angles)))
    counts_pro = count_torsion_angles(phi_psi_angles=phi_psi_angles,
                                      resolution=resolution,
                                      normalized=True)
    # pre-PRO
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_prepro()
    print(
        "Couting torsions for pre-PRO, {} samples, this might take a while...".
        format(len(phi_psi_angles)))
    counts_prepro = count_torsion_angles(phi_psi_angles=phi_psi_angles,
                                         resolution=resolution,
                                         normalized=True)
    # General
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_general(
    )
    print(
        "Couting torsions for General, {} samples, this might take a while...".
        format(len(phi_psi_angles)))
    counts_general = count_torsion_angles(phi_psi_angles=phi_psi_angles,
                                          resolution=resolution,
                                          normalized=True)

    if save_file_path is not None:

        dir_path = os.path.dirname(save_file_path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        np.savez(save_file_path,
                 gly=counts_gly,
                 pro=counts_pro,
                 prepro=counts_prepro,
                 general=counts_general)

    return (counts_gly, counts_pro, counts_prepro, counts_general)


def compute_gaussian_kde_densities_from_directory(
    dir_path: str,
    save_file_path: Optional[str] = None,
    b_factor_threshold: float = 30,
    resolution: int = 180,
    num_processes: int = 4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    residue_torsion_collection = read_residue_torsion_collection_from_directory(
        dir_path=dir_path,
        b_factor_threshold=b_factor_threshold,
        num_processes=num_processes)

    # Computing Gaussian KDE
    # This is a single-thread implementation and is somewhat slow.
    # GLY
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_gly()
    print(
        "Computing Gaussian kernel for GLY, {} samples, this might take a while..."
        .format(len(phi_psi_angles)))
    gaussian_density_gly = compute_gaussian_kde_density(
        phi_psi_angles=phi_psi_angles, resolution=resolution)
    # PRO
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_pro()
    print(
        "Computing Gaussian kernel for PRO, {} samples, this might take a while..."
        .format(len(phi_psi_angles)))
    gaussian_density_pro = compute_gaussian_kde_density(
        phi_psi_angles=phi_psi_angles, resolution=resolution)
    # pre-PRO
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_prepro()
    print(
        "Computing Gaussian kernel for pre-PRO, {} samples, this might take a while..."
        .format(len(phi_psi_angles)))
    gaussian_density_prepro = compute_gaussian_kde_density(
        phi_psi_angles=phi_psi_angles, resolution=resolution)
    # General
    phi_psi_angles = residue_torsion_collection.collect_torsion_angles_general(
    )
    print(
        "Computing Gaussian kernel for general, {} samples, this might take a while..."
        .format(len(phi_psi_angles)))
    gaussian_density_general = compute_gaussian_kde_density(
        phi_psi_angles=phi_psi_angles, resolution=resolution)

    # print(np.sum(gaussian_density_gly))
    # print(np.sum(gaussian_density_pro))
    # print(np.sum(gaussian_density_prepro))
    # print(np.sum(gaussian_density_general))

    if save_file_path is not None:

        dir_path = os.path.dirname(save_file_path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        np.savez(save_file_path,
                 gly=gaussian_density_gly,
                 pro=gaussian_density_pro,
                 prepro=gaussian_density_prepro,
                 general=gaussian_density_general)

    return (gaussian_density_gly, gaussian_density_pro,
            gaussian_density_prepro, gaussian_density_general)
