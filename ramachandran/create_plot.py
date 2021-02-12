from typing import List, Tuple, Optional
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from ramachandran.read_pdb import read_protein
from ramachandran.compute_dihedral_angle import protein_backbone_dihedral_angle_phi, protein_backbone_dihedral_angle_psi


def collect_dihedral_angles(pdb_filepath: str) -> List[Tuple[float, float]]:

    protein = read_protein(pdb_filepath=pdb_filepath)

    dihedrals = []

    for chain_identifier in protein:
        chain = protein[chain_identifier]
        for residue_sequence_number in chain:
            residue = chain[residue_sequence_number]
            # Skip the first, the last, and problematic residues
            if residue_sequence_number - 1 not in chain or residue_sequence_number + 1 not in chain:
                continue
            last_residue = chain[residue_sequence_number - 1]
            next_residue = chain[residue_sequence_number + 1]
            n = residue["N"]
            c_alpha = residue["CA"]
            c = residue["C"]
            c_minus = last_residue["C"]
            n_plus = next_residue["N"]

            phi = protein_backbone_dihedral_angle_phi(c_minus=c_minus, n=n, c_alpha=c_alpha, c=c)
            psi = protein_backbone_dihedral_angle_psi(n=n, c_alpha=c_alpha, c=c, n_plus=n_plus)

            phi = np.rad2deg(phi)
            psi = np.rad2deg(psi)

            dihedrals.append((phi, psi))
    
    return dihedrals


def create_ramachandran_plot(pdb_filepath: str, plot_filepath: str, protein_name: Optional[str] = None) -> None:

    dihedrals = collect_dihedral_angles(pdb_filepath=pdb_filepath)

    x = []
    y = []

    for phi, psi in dihedrals:
        x.append(phi)
        y.append(psi)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=10)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))

    ax.grid(linestyle='--')


    ax.set_xlabel(r'${\phi}$', fontsize=16, fontweight='bold')
    ax.set_ylabel(r'${\psi}$', fontsize=16, fontweight='bold')


    # ax.legend()
    fig.savefig(plot_filepath, format='svg', dpi=600, bbox_inches='tight')
    plt.close()
            

    return