from typing import Iterable

from rdkit import Chem


def calculate_ideal_internal_angle_sum(ring_size: int) -> float:
    """Return the ideal internal angle sum for a ring of given size."""
    return (ring_size - 2) * 180 / ring_size


def calculate_strain_factor(
    mol: Chem.Mol,
    rings: Iterable[Iterable[int]],
    lambda_macrocycle: float = 0.5,
    macrocycle_threshold: int = 12,
) -> float:
    """Calculate strain factor with reduced contribution for macrocycles."""
    strain_factor = 0.0
    for ring in rings:
        ring_size = len(ring)
        ideal_angle_sum = calculate_ideal_internal_angle_sum(ring_size)
        lambda_M = lambda_macrocycle if ring_size >= macrocycle_threshold else 1.0
        strain_factor += (360 / (360 - ideal_angle_sum)) / ring_size * lambda_M
    return strain_factor


def count_macrocycles(mol: Chem.Mol, min_ring_size: int = 12) -> int:
    """Return number of rings with size >= ``min_ring_size``."""
    if mol is None:
        return 0
    sssr = Chem.GetSymmSSSR(mol)
    count = 0
    for ring in sssr:
        if len(ring) >= min_ring_size:
            count += 1
    return count
