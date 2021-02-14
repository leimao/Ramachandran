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
    sequence = {}

    with open(pdb_filepath, "r") as fp:

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
                    raise RuntimeError("Atom serial number has to be integer. Got {}".format(items[1]))

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
                
                if "is_pre_proline" not in protein[chain_identifier][
                        residue_sequence_number]:
                    
                    chain_sequence = sequence[chain_identifier]
                    num_residues = len(chain_sequence)

                    if residue_sequence_number == num_residues:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = False
                    else:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = (chain_sequence[residue_sequence_number] == "PRO")

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

    return protein


def read_pdbx(pdbx_filepath: str) -> Dict[Any, Any]:

    sequence = {}

    with open(pdbx_filepath, "r") as fp:

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
    with open(pdbx_filepath, "r") as fp:

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
                    raise RuntimeError("Atom serial number has to be integer. Got {}".format(pdbx_filepath))

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

                if "is_pre_proline" not in protein[chain_identifier][
                        residue_sequence_number]:
                    
                    chain_sequence = sequence[chain_sequence_identifier]
                    num_residues = len(chain_sequence)

                    if residue_sequence_number == num_residues:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = False
                    else:
                        protein[chain_identifier][residue_sequence_number][
                            "is_pre_proline"] = (chain_sequence[residue_sequence_number] == "PRO")

                protein[chain_identifier][residue_sequence_number][
                    atom_name] = atom_coodinates

            if line.strip() == "_atom_site.pdbx_PDB_model_num":
                is_atom = True

    # print(protein)
    return protein
