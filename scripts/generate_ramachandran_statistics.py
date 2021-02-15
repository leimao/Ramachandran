import ramachandran.statistics

ramachandran.statistics.count_torsion_angles_from_directory(dir_path="./pdbx_collections_resolution_1p0", save_file_path="./data/probability.npz", b_factor_threshold=30, resolution=90, num_processes=12)

ramachandran.statistics.compute_gaussian_kde_densities_from_directory(dir_path="./pdbx_collections_resolution_1p0", save_file_path="./data/gaussian_density.npz", b_factor_threshold=30, resolution=360, num_processes=12)