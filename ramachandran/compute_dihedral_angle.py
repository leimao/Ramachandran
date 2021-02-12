from collections.abc import Sequence
import numpy as np
from numpy import linalg as LA


def polymer_physics_dihedral_angle(a: Sequence, b: Sequence, c: Sequence) -> float:
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
    sin_phi = np.inner(b, np.cross(a=axb, b=bxc)) / (norm_b * norm_axb * norm_bxc)

    angle = np.arctan2(sin_phi, cos_phi)

    return angle


def protein_backbone_dihedral_angle_phi(c_minus: Sequence, n: Sequence, c_alpha: Sequence, c: Sequence) -> float:
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


def protein_backbone_dihedral_angle_psi(n: Sequence, c_alpha: Sequence, c: Sequence, n_plus) -> float:
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
