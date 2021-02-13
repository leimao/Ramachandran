from typing import List, Tuple, Optional
import os
import math
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
import numpy as np
from scipy import stats
from ramachandran.read_structure import read_pdb, read_pdbx
from ramachandran.compute_dihedral_angle import protein_backbone_dihedral_angle_phi, protein_backbone_dihedral_angle_psi, collect_dihedral_angles


def count_dihedral_angles(dihedral_angles: List[Tuple[float, float]],
                          num_bins_per_dimension: int = 180) -> np.ndarray:

    counts = np.zeros(shape=(num_bins_per_dimension, num_bins_per_dimension),
                      dtype=np.uint32)

    for dihedral_angle in dihedral_angles:

        phi, psi = dihedral_angle

        # i = int((phi + 180) / 360 * num_bins_per_dimension)
        # j = int((psi + 180) / 360 * num_bins_per_dimension)
        i = int((180 - psi) / 360 * num_bins_per_dimension)
        j = int((phi + 180) / 360 * num_bins_per_dimension)

        # If i or j = num_bins_per_dimension, adjust it.
        if i == num_bins_per_dimension:
            i = num_bins_per_dimension - 1
        if j == num_bins_per_dimension:
            j = num_bins_per_dimension - 1

        counts[i, j] += 1

    return counts


def count_dihedral_angles_from_file(filepath: str,
                                    num_bins_per_dimension: int = 180,
                                    b_factor_threshold: float = 30
                                    ) -> np.ndarray:

    try:
        dihedral_angles = collect_dihedral_angles(filepath=filepath, b_factor_threshold=b_factor_threshold)
    except RuntimeError:
        dihedral_angles = []

    counts = count_dihedral_angles(
        dihedral_angles=dihedral_angles,
        num_bins_per_dimension=num_bins_per_dimension)

    return counts


def count_dihedral_angles_from_files(filepaths: List[str],
                                     num_bins_per_dimension: int = 180,
                                     num_processes: int = 4,
                                     b_factor_threshold: float = 30) -> np.ndarray:

    pool = Pool(processes=num_processes)

    partial_func = partial(count_dihedral_angles_from_file,
                           num_bins_per_dimension=num_bins_per_dimension,
                           b_factor_threshold=b_factor_threshold)

    jobs = [
        pool.apply_async(func=partial_func, args=(filepath, ))
        for filepath in filepaths
    ]

    pool.close()

    results = []
    for job in tqdm(jobs):
        results.append(job.get())

    # Reduce the counts
    reduced_counts = sum(results)

    return reduced_counts


def count_dihedral_angles_from_directory(dir_path: str,
                                         data_filepath: Optional[str] = None,
                                         max_batch_size: int = 1000,
                                         num_bins_per_dimension: int = 180,
                                         b_factor_threshold: float = 30,
                                         num_processes: int = 4) -> np.ndarray:

    valid_structure_filepaths = []

    for f in os.listdir(dir_path):

        _, file_extension = os.path.splitext(f)
        filepath = os.path.join(dir_path, f)

        if os.path.isfile(filepath) and (file_extension == ".pdb"
                                         or file_extension == ".cif"):
            valid_structure_filepaths.append(filepath)

    # To limit the usage of memory, we process the files in batches.
    num_batches = int(
        math.ceil(len(valid_structure_filepaths) / max_batch_size))

    reduced_counts = np.zeros(shape=(num_bins_per_dimension,
                                     num_bins_per_dimension),
                              dtype=np.uint32)

    for i in range(num_batches):

        print("Batch: {:03d}/{:03d}".format(i + 1, num_batches))

        filepaths = valid_structure_filepaths[i * max_batch_size:(i + 1) *
                                              max_batch_size]

        new_reduced_counts = count_dihedral_angles_from_files(
            filepaths=filepaths,
            num_bins_per_dimension=num_bins_per_dimension,
            num_processes=num_processes,
            b_factor_threshold=b_factor_threshold)

        reduced_counts = sum([reduced_counts, new_reduced_counts])

    counts_density = reduced_counts.astype(np.float32) / np.sum(reduced_counts)

    if data_filepath is not None:

        np.savez(data_filepath, density_map=counts_density)

    return counts_density



def collect_dihedral_angles_from_file(filepath: str,
                                    b_factor_threshold: float = 30
                                    ) -> List[Tuple[float, float]]:

    try:
        dihedral_angles = collect_dihedral_angles(filepath=filepath, b_factor_threshold=b_factor_threshold)
    except RuntimeError:
        dihedral_angles = []

    return dihedral_angles


def collect_dihedral_angles_from_files(filepaths: List[str],
                                     num_processes: int = 4,
                                     b_factor_threshold: float = 30) -> List[Tuple[float, float]]:

    pool = Pool(processes=num_processes)

    partial_func = partial(collect_dihedral_angles_from_file,
                           b_factor_threshold=b_factor_threshold)

    jobs = [
        pool.apply_async(func=partial_func, args=(filepath, ))
        for filepath in filepaths
    ]

    pool.close()

    dihedral_angles = []
    for job in tqdm(jobs):
        dihedral_angles += job.get()

    return dihedral_angles

def collect_dihedral_angles_from_directory(dir_path: str,
                                         b_factor_threshold: float = 30,
                                         num_processes: int = 4) -> List[Tuple[float, float]]:

    valid_structure_filepaths = []

    for f in os.listdir(dir_path):

        _, file_extension = os.path.splitext(f)
        filepath = os.path.join(dir_path, f)

        if os.path.isfile(filepath) and (file_extension == ".pdb"
                                         or file_extension == ".cif"):
            valid_structure_filepaths.append(filepath)

    dihedral_angles = collect_dihedral_angles_from_files(filepaths=valid_structure_filepaths, num_processes=num_processes, b_factor_threshold=b_factor_threshold)

    return dihedral_angles

def compute_gaussian_kde_density(dihedral_angles: List[Tuple[float, float]], resolution: int = 180) -> np.ndarray:

    # Turn the list of tuples to numpy array
    # Shape: [2, len(dihedral_angles)]
    values = np.array(dihedral_angles).T

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

def compute_gaussian_kde_density_from_directory(dir_path: str,
                                         data_filepath: Optional[str] = None,
                                         b_factor_threshold: float = 30,
                                         resolution: int = 180,
                                         num_processes: int = 4) -> np.ndarray:


    dihedral_angles = collect_dihedral_angles_from_directory(dir_path=dir_path,
                                         b_factor_threshold=b_factor_threshold,
                                         num_processes=num_processes)

    # Computing Gaussian KDE
    # This is a single-thread implementation and is somewhat slow.
    print("Computing Gaussian kernel, this might take a while...")
    gaussian_density = compute_gaussian_kde_density(dihedral_angles=dihedral_angles, resolution=resolution)

    if data_filepath is not None:

        np.savez(data_filepath, density_map=gaussian_density)

    return gaussian_density