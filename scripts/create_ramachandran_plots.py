#!/usr/bin/python3

import argparse
import ramachandran.plot


if __name__ == "__main__":

    ramachandran.plot.create_ramachandran_plots_from_file(file_path="./2hik.pdb",
                                save_dir_path="./2hik_pdb",
                                protein_name="2HIK")

    ramachandran.plot.create_ramachandran_plots_from_file(file_path="./2hik.cif",
                                save_dir_path="./2hik_cif",
                                protein_name="2HIK")