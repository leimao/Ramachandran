import ramachandran.read_pdb
import ramachandran.create_plot

protein = ramachandran.read_pdb.read_protein(pdb_filepath="./2hik.pdb")

ramachandran.create_plot.create_ramachandran_plot(pdb_filepath="./2hik.pdb", plot_filepath="./2hik.svg")