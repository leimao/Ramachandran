import ramachandran.compute_trust_regions
import numpy as np

# counts = ramachandran.compute_trust_regions.count_dihedral_angles_from_directory(dir_path="./pdbx_collections_resolution_1p2", data_filepath="./general_region.npz", num_bins_per_dimension=180, b_factor_threshold=30, num_processes=8)

# print(counts)

# dihedral_angles = ramachandran.compute_trust_regions.collect_dihedral_angles_from_directory(dir_path="./pdbx_collections_resolution_1p0", data_filepath="./general_region.npz", b_factor_threshold=30, num_processes=8)

# print(dihedral_angles)


# print(np.array(dihedral_angles))

density = ramachandran.compute_trust_regions.compute_gaussian_kde_density_from_directory(dir_path="./pdbx_collections_resolution_1p0",
                                         data_filepath="./general_region_gaussian.npz",
                                         b_factor_threshold=30,
                                         resolution=360,
                                         num_processes=12)

print(density.shape)