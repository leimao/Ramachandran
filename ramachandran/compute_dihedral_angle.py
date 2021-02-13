from collections.abc import Sequence
from typing import List, Tuple, Optional
import os
import numpy as np
from numpy import linalg as LA
from ramachandran.read_structure import read_pdb, read_pdbx


def polymer_physics_dihedral_angle(a: Sequence, b: Sequence,
                                   c: Sequence) -> float:
    """Given a plane defined by vectors a and b, and a plane defined by vectors b and c, compute the dihedral angle of the two planes. The dihedral angle is in a range of [-pi, pi].

    Args:
        a (Sequence): [description]
        b (Sequence): [description]
        c (Sequence): [description]

    Returns:
        float: [description]
    """

    # https://en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics

    axb = np.cross(a=a, b=b)
    norm_axb = LA.norm(axb)
    bxc = np.cross(a=b, b=c)
    norm_bxc = LA.norm(bxc)

    norm_b = LA.norm(b)

    cos_phi = np.inner(axb, bxc) / (norm_axb * norm_bxc)

    # We use the definition of cross product to compute sin phi.
    # axb, b, bxc is a right handed system.
    sin_phi = np.inner(b, np.cross(a=axb,
                                   b=bxc)) / (norm_b * norm_axb * norm_bxc)

    angle = np.arctan2(sin_phi, cos_phi)

    return angle


def protein_backbone_dihedral_angle_phi(c_minus: Sequence, n: Sequence,
                                        c_alpha: Sequence,
                                        c: Sequence) -> float:
    """Given the atom coordinates of C_(-1), N, C_alpha, and C, compute the dihedral angle phi.

    Args:
        c_minus (Sequence): [description]
        n (Sequence): [description]
        c_alpha (Sequence): [description]
        c (Sequence): [description]

    Returns:
        float: [description]
    """

    # Compute three vectors
    b1 = np.array(n) - np.array(c_minus)
    b2 = np.array(c_alpha) - np.array(n)
    b3 = np.array(c) - np.array(c_alpha)

    angle = polymer_physics_dihedral_angle(a=b1, b=b2, c=b3)

    return angle


def protein_backbone_dihedral_angle_psi(n: Sequence, c_alpha: Sequence,
                                        c: Sequence, n_plus) -> float:
    """Given the atom coordinates of N, C_alpha, C, and N_(+1) compute the dihedral angle psi.

    Args:
        n (Sequence): [description]
        c_alpha (Sequence): [description]
        c (Sequence): [description]
        n_plus ([type]): [description]

    Returns:
        float: [description]
    """

    # Compute three vectors
    b1 = np.array(c_alpha) - np.array(n)
    b2 = np.array(c) - np.array(c_alpha)
    b3 = np.array(n_plus) - np.array(c)

    angle = polymer_physics_dihedral_angle(a=b1, b=b2, c=b3)

    return angle


def collect_dihedral_angles(filepath: str, b_factor_threshold: Optional[float] = None) -> List[Tuple[float, float]]:

    _, file_extension = os.path.splitext(filepath)

    if file_extension == ".pdb":
        protein = read_pdb(pdb_filepath=filepath)
    elif file_extension == ".cif":
        protein = read_pdbx(pdbx_filepath=filepath)
    else:
        raise RuntimeError(
            "Only files with extensions of pdb and cif are supported.")

    dihedral_angles = []

    for chain_identifier in protein:
        chain = protein[chain_identifier]
        for residue_sequence_number in chain:
            residue = chain[residue_sequence_number]
            # Skip the first, the last, and problematic residues
            if residue_sequence_number - 1 not in chain or residue_sequence_number + 1 not in chain:
                continue
            last_residue = chain[residue_sequence_number - 1]
            next_residue = chain[residue_sequence_number + 1]

            # Skip the residues that has missing information to compute dihedral angles
            if "N" not in residue or "CA" not in residue or "C" not in residue or "C" not in last_residue or "N" not in next_residue:
                continue

            n, b_factor_n = residue["N"]
            c_alpha, b_factor_c_alpha = residue["CA"]
            c, b_factor_c = residue["C"]
            c_minus, b_factor_c_minus = last_residue["C"]
            n_plus, b_factor_n_plus = next_residue["N"]

            if b_factor_threshold is not None:
                if b_factor_n > b_factor_threshold or b_factor_c_alpha > b_factor_threshold or b_factor_c > b_factor_threshold or b_factor_c_minus > b_factor_threshold or b_factor_n_plus > b_factor_threshold:
                    continue

            phi = protein_backbone_dihedral_angle_phi(c_minus=c_minus,
                                                      n=n,
                                                      c_alpha=c_alpha,
                                                      c=c)
            psi = protein_backbone_dihedral_angle_psi(n=n,
                                                      c_alpha=c_alpha,
                                                      c=c,
                                                      n_plus=n_plus)

            phi = np.rad2deg(phi)
            psi = np.rad2deg(psi)

            dihedral_angles.append((phi, psi))

    return dihedral_angles
