from typing import Any, Dict


def read_pdb(pdb_filepath: str) -> Dict[Any, Any]:
    """Read the basic information of a protein from a PDB file. 
        The function does minimal format checking for the PDB format.
        Please make sure the PDB format is correct before calling the function.

    Args:
        pdb_filepath (str): [description]

    Raises:
        RuntimeError: [description]
        RuntimeError: [description]
        RuntimeError: [description]

    Returns:
        Dict[Any, Any]: [description]
    """

    protein = {}

    with open(pdb_filepath, "r") as fp:
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
                    raise RuntimeError("Atom serial number has to be integer.")

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
                    raise RuntimeError(
                        "B factor has to be a real value.")

                atom_coodinates = ((x, y, z), b_factor)

                if chain_identifier not in protein:
                    protein[chain_identifier] = {}

                if residue_sequence_number not in protein[chain_identifier]:
                    protein[chain_identifier][residue_sequence_number] = {}

                if "residue" not in protein[chain_identifier][
                        residue_sequence_number]:
                    protein[chain_identifier][residue_sequence_number][
                        "residue"] = residue_name

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

    return protein


def read_pdbx(pdbx_filepath: str) -> Dict[Any, Any]:

    protein = {}

    # http://ww1.iucr.org/iucr-top/cif/mmcif/workshop/mmCIF-tutorials/intro/atom.htm
    with open(pdbx_filepath, "r") as fp:
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
                    raise RuntimeError("Atom serial number has to be integer.")

                atom_name = items[3]
                residue_name = items[5]
                chain_identifier = items[6]
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
                    raise RuntimeError(
                        "B factor has to be a real value.")

                atom_coodinates = ((x, y, z), b_factor)

                if chain_identifier not in protein:
                    protein[chain_identifier] = {}

                if residue_sequence_number not in protein[chain_identifier]:
                    protein[chain_identifier][residue_sequence_number] = {}

                if "residue" not in protein[chain_identifier][
                        residue_sequence_number]:
                    protein[chain_identifier][residue_sequence_number][
                        "residue"] = residue_name

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

    return protein
