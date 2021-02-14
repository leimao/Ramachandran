from typing import Any, Dict, Optional, List
import os
import numpy as np
from ramachandran.geometry import protein_backbone_dihedral_angle_phi, protein_backbone_dihedral_angle_psi
from ramachandran.torsion import ResidueTorsionCollection, ResidueTorsion
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial


def read_pdb(pdb_file_path: str) -> Dict[Any, Any]:
    """Read the basic information of a protein from a PDB file. 
        The function does minimal format checking for the PDB format.
        Please make sure the PDB format is correct before calling the function.

    Args:
        pdb_file_path (str): file path to the PDB file.

    Returns:
        Dict[Any, Any]: A dictionary containing the basic protein information.
    """
    sequence = {}

    with open(pdb_file_path, "r") as fp:

        for line in fp:
            # Check if the line is SEQRES
            if line[0:6] == "SEQRES":
                items = line.split()

                chain_identifier = items[2]
                amino_acids = items[4:]

                if chain_identifier not in sequence:
                    sequence[chain_identifier] = []
                sequence[chain_identifier] += amino_acids

    protein = {}

    with open(pdb_file_path, "r") as fp:

        for line in fp:
            # Check if the line is ATOM
            # 1-4 Record Type
            if line[0:4] == "ATOM":
                # It is not best practice for dealing with PDB format.
                # Use it for the ATOM section for now.
                items = line.split()

                try:
                    atom_serial_number = int(items[1])
                except:
                    raise RuntimeError(
                        "Atom serial number has to be integer. Got {}".format(
                            items[1]))

                atom_name = items[2]
                residue_name = items[3]
                chain_identifier = items[4]
                try:
                    residue_sequence_number = int(items[5])
                except:
                    raise RuntimeError(
                        "Residue sequence number has to be integer.")

                try:
                    x = float(items[6])
                    y = float(items[7])
                    z = float(items[8])
                except:
                    raise RuntimeError(
                        "Atom coordinate has to be a real value.")

                try:
                    b_factor = float(items[10])
                except:
                    raise RuntimeError("B factor has to be a real value.")

                atom_coodinates = ((x, y, z), b_factor)

                if chain_identifier not in protein:
                    protein[chain_identifier] = {}

                if residue_sequence_number not in protein[chain_identifier]:
                    protein[chain_identifier][residue_sequence_number] = {}

                if "residue" not in protein[chain_identifier][
                        residue_sequence_number]:
                    protein[chain_identifier][residue_sequence_number][
                        "residue"] = residue_name

                if "is_pre_proline" not in protein[chain_identifier][
                        residue_sequence_number]:

                    chain_sequence = sequence[chain_identifier]
                    num_residues = len(chain_sequence)

                    if residue_sequence_number == num_residues:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = False
                    else:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = (
                                chain_sequence[residue_sequence_number] ==
                                "PRO")

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

    return protein


def read_pdbx(pdbx_file_path: str) -> Dict[Any, Any]:
    """Read the basic information of a protein from a PDBx file. 
        The function does minimal format checking for the PDBx format.
        Please make sure the PDBx format is correct before calling the function.

    Args:
        pdbx_file_path (str): file path to the PDBx file.

    Returns:
        Dict[Any, Any]: A dictionary containing the basic protein information.
    """

    sequence = {}

    with open(pdbx_file_path, "r") as fp:

        is_sequence = False

        for line in fp:

            if line.strip() == "#":
                is_sequence = False

            # Check if the line is SEQRES
            if is_sequence:
                items = line.split()
                assert len(items) == 4, line

                chain_identifier = items[0]
                amino_acid = items[2]

                if chain_identifier not in sequence:
                    sequence[chain_identifier] = []
                sequence[chain_identifier].append(amino_acid)

            if line.strip() == "_entity_poly_seq.hetero":
                is_sequence = True

    protein = {}

    # http://ww1.iucr.org/iucr-top/cif/mmcif/workshop/mmCIF-tutorials/intro/atom.htm
    with open(pdbx_file_path, "r") as fp:

        is_atom = False

        for line in fp:

            if line.strip() == "#":
                is_atom = False

            # Check if the line is ATOM
            # 1-4 Record Type
            if is_atom and line[0:4] == "ATOM":
                # It is not best practice for dealing with PDB format.
                # Use it for the ATOM section for now.
                items = line.split()

                try:
                    atom_serial_number = int(items[1])
                except:
                    raise RuntimeError(
                        "Atom serial number has to be integer. Got {}".format(
                            pdbx_file_path))

                atom_name = items[3]
                residue_name = items[5]
                chain_identifier = items[6]
                chain_sequence_identifier = items[7]

                try:
                    residue_sequence_number = int(items[8])
                except:
                    raise RuntimeError(
                        "Residue sequence number has to be integer.")

                try:
                    x = float(items[10])
                    y = float(items[11])
                    z = float(items[12])
                except:
                    raise RuntimeError(
                        "Atom coordinate has to be a real value.")

                try:
                    b_factor = float(items[14])
                except:
                    raise RuntimeError("B factor has to be a real value.")

                atom_coodinates = ((x, y, z), b_factor)

                if chain_identifier not in protein:
                    protein[chain_identifier] = {}

                if residue_sequence_number not in protein[chain_identifier]:
                    protein[chain_identifier][residue_sequence_number] = {}

                if "residue" not in protein[chain_identifier][
                        residue_sequence_number]:
                    protein[chain_identifier][residue_sequence_number][
                        "residue"] = residue_name

                if "is_pre_proline" not in protein[chain_identifier][
                        residue_sequence_number]:

                    chain_sequence = sequence[chain_sequence_identifier]
                    num_residues = len(chain_sequence)

                    if residue_sequence_number == num_residues:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = False
                    else:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = (
                                chain_sequence[residue_sequence_number] ==
                                "PRO")

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

            if line.strip() == "_atom_site.pdbx_PDB_model_num":
                is_atom = True

    return protein


def read_residue_torsion_collection_from_protein(
        protein: Dict[Any, Any],
        b_factor_threshold: Optional[float] = None
) -> ResidueTorsionCollection:
    """Read torsion collection from a protein dictionary.

    Args:
        protein (Dict[Any, Any]): A protein dictionary.
        b_factor_threshold (Optional[float], optional): B-factor indicating the uncertainty of the atom coordinates. Defaults to None.

    Returns:
        ResidueTorsionCollection: A residue torsion collection for the protein given.
    """

    residue_torsion_collection = ResidueTorsionCollection()

    for chain_identifier in protein:

        chain = protein[chain_identifier]
        for residue_sequence_number in chain:
            residue = chain[residue_sequence_number]
            residue_name = residue["residue"]
            # Skip the first, the last, and problematic residues
            if residue_sequence_number - 1 not in chain or residue_sequence_number + 1 not in chain:
                continue
            last_residue = chain[residue_sequence_number - 1]
            next_residue = chain[residue_sequence_number + 1]
            next_residue_name = next_residue["residue"]

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
            is_pre_proline = (next_residue_name == "PRO")

            residue_torsion = ResidueTorsion(phi=phi,
                                             psi=psi,
                                             residue_type=residue_name,
                                             is_pre_proline=is_pre_proline)

            residue_torsion_collection.append(residue_torsion=residue_torsion)

    return residue_torsion_collection


def read_residue_torsion_collection_from_file(
        file_path: str,
        b_factor_threshold: Optional[float] = None
) -> ResidueTorsionCollection:
    """Read torsion collection from a PDB or PDBx file.

    Args:
        file_path (str): file path to a PDB or PDBx file.
        b_factor_threshold (Optional[float], optional): B-factor indicating the uncertainty of the atom coordinates. Defaults to None.

    Raises:
        RuntimeError: When non-PDB or non-PDBx files were provided, raise RuntimeError.

    Returns:
        ResidueTorsionCollection: A residue torsion collection for the PDB or PDBx file given.
    """

    _, file_extension = os.path.splitext(file_path)

    if file_extension == ".pdb":
        protein = read_pdb(pdb_file_path=file_path)
    elif file_extension == ".cif":
        protein = read_pdbx(pdbx_file_path=file_path)
    else:
        raise RuntimeError(
            "Only files with extensions of pdb and cif are supported.")

    residue_torsion_collection = read_residue_torsion_collection_from_protein(
        protein=protein, b_factor_threshold=b_factor_threshold)

    return residue_torsion_collection


def read_residue_torsion_collection_from_file_resolved(
        file_path: str,
        b_factor_threshold: Optional[float] = None
) -> ResidueTorsionCollection:
    try:
        residue_torsion_collection = read_residue_torsion_collection_from_file(
            file_path=file_path, b_factor_threshold=b_factor_threshold)
    except:
        residue_torsion_collection = ResidueTorsionCollection()
    return residue_torsion_collection


def read_residue_torsion_collection_from_files(
        file_paths: List[str],
        num_processes: int = 4,
        b_factor_threshold: Optional[float] = None
) -> ResidueTorsionCollection:

    pool = Pool(processes=num_processes)

    residue_torsion_collection_all = ResidueTorsionCollection()

    partial_func = partial(read_residue_torsion_collection_from_file_resolved,
                           b_factor_threshold=b_factor_threshold)

    jobs = [
        pool.apply_async(func=partial_func, args=(file_path, ))
        for file_path in file_paths
    ]

    pool.close()

    for job in tqdm(jobs):
        residue_torsion_collection_all += job.get()

    return residue_torsion_collection_all


def read_residue_torsion_collection_from_directory(
        dir_path: str,
        b_factor_threshold: float = 30,
        num_processes: int = 4) -> ResidueTorsionCollection:

    valid_file_paths = []

    for f in os.listdir(dir_path):

        _, file_extension = os.path.splitext(f)
        file_path = os.path.join(dir_path, f)

        if os.path.isfile(file_path) and (file_extension == ".pdb"
                                          or file_extension == ".cif"):
            valid_file_paths.append(file_path)

    residue_torsion_collection = read_residue_torsion_collection_from_files(
        file_paths=valid_file_paths,
        num_processes=num_processes,
        b_factor_threshold=b_factor_threshold)

    return residue_torsion_collection
