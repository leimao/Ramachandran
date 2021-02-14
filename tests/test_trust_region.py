import ramachandran.compute_trust_regions
import numpy as np

ramachandran.compute_trust_regions.compute_gaussian_kde_densities_from_directory(dir_path="./pdbx_collections_resolution_1p0",
                                         data_filepath="./general_region_gaussian.npz",
                                         b_factor_threshold=30,
                                         resolution=180,
                                         num_processes=12)

