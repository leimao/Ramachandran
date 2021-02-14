from __future__ import annotations
from typing import List, Optional, Sequence, Set, Tuple

RESIDUE_TYPES = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
}


class ResidueTorsion(object):
    """Data structure for storing residue torsion information.
    """
    def __init__(self,
                 phi: float,
                 psi: float,
                 omega: float = 180,
                 residue_type: Optional[str] = None,
                 is_pre_proline: Optional[bool] = None) -> None:

        self.phi = phi
        self.psi = psi
        self.omega = omega
        self.residue_type = residue_type
        self.is_pre_proline = is_pre_proline

    def _sanity_check(self):

        if self.phi < -180 or self.phi > 180:
            raise ValueError(
                f"Torsion angle phi is defined in a range of [-180, 180]. Got {self.phi}"
            )

        if self.psi < -180 or self.psi > 180:
            raise ValueError(
                f"Torsion angle psi is defined in a range of [-180, 180]. Got {self.psi}"
            )

        if self.omega < -180 or self.omega > 180:
            raise ValueError(
                f"Torsion angle omega is defined in a range of [-180, 180]. Got {self.omega}"
            )

        if self.residue_type is not None and self.residue_type not in RESIDUE_TYPES:
            raise ValueError(
                f"Unsupported amino acid type. Got {self.residue_type}")


class ResidueTorsionCollection(object):
    """A collection of ResidueTorsion objects for analysis.
    """
    def __init__(
            self,
            residue_torsions: Optional[List[ResidueTorsion]] = None) -> None:

        # The ResidueTorsionCollection object will own a copy of the data.
        if residue_torsions is None:
            self.residue_torsions = []
        else:
            self.residue_torsions = residue_torsions.copy()

    def append(self, residue_torsion: ResidueTorsion) -> None:

        self.residue_torsions.append(residue_torsion)

    def collect_torsion_angles_by_residue_types(
            self,
            residue_types: Set[str],
            is_pre_proline: Optional[bool] = None
    ) -> List[Tuple[float, float]]:
        """Collect torsion angles phi and psi by residue type.

        Args:
            residue_type (Set[str]): A sequence of residue types that needs to be collected.

        Returns:
            List[Tuple[float, float]]: A list of phi psi angle tuples.
        """

        phi_psi_angles = []

        for residue_torsion in self.residue_torsions:
            if residue_torsion.residue_type in residue_types:
                if is_pre_proline is not None and residue_torsion.is_pre_proline != is_pre_proline:
                    continue
                phi_psi_angles.append(
                    (residue_torsion.phi, residue_torsion.psi))

        return phi_psi_angles

    def collect_torsion_angles_pro(self) -> List[Tuple[float, float]]:

        return self.collect_torsion_angles_by_residue_types(
            residue_types={"PRO"}, is_pre_proline=False)

    def collect_torsion_angles_prepro(self) -> List[Tuple[float, float]]:

        return self.collect_torsion_angles_by_residue_types(
            residue_types=RESIDUE_TYPES, is_pre_proline=True)

    def collect_torsion_angles_gly(self) -> List[Tuple[float, float]]:

        return self.collect_torsion_angles_by_residue_types(
            residue_types={"GLY"}, is_pre_proline=False)

    def collect_torsion_angles_general(self) -> List[Tuple[float, float]]:

        general_residue_types = RESIDUE_TYPES.copy()
        general_residue_types.remove("PRO")
        general_residue_types.remove("GLY")

        return self.collect_torsion_angles_by_residue_types(
            residue_types=general_residue_types, is_pre_proline=False)

    # Override +=
    def __iadd__(self,
                 other: ResidueTorsionCollection) -> ResidueTorsionCollection:

        self.residue_torsions.extend(other.residue_torsions)

        return self

    # Override +
    def __add__(self,
                other: ResidueTorsionCollection) -> ResidueTorsionCollection:

        return ResidueTorsionCollection(
            residue_torsions=self.residue_torsions + other.residue_torsions)
