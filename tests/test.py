import ramachandran.read_structure
import ramachandran.create_plot

protein = ramachandran.read_structure.read_pdb(pdb_filepath="./2hik.pdb")

ramachandran.create_plot.create_ramachandran_plot(filepath="./2hik.cif", plot_filepath="./2hik_cif.svg")

ramachandran.create_plot.create_ramachandran_plot(filepath="./2hik.pdb", plot_filepath="./2hik.svg")

#ramachandran.create_plot.create_ramachandran_plot(filepath="./2hik.pdb", plot_filepath="./2hik.svg")

ramachandran.create_plot.create_ramachandran_trusted_regsion_plot(trusted_regsion_filepath="./general_region_gaussian.npz", plot_filepath="./test.svg")