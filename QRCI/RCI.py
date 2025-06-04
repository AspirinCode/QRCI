from typing import Union
from rdkit import Chem


def calculate_RCI(smiles_or_mol: Union[str, Chem.Mol]) -> float | None:
    """Calculate the Ring Complexity Index (RCI).

    Parameters
    ----------
    smiles_or_mol : str | Chem.Mol
        Molecule given as a SMILES string or RDKit Mol object.

    Returns
    -------
    float | None
        The RCI value or ``None`` if the molecule could not be parsed.
    """
    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
        if mol is None:
            return None
    else:
        mol = smiles_or_mol

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return 0.0

    TRS = sum(len(ring) for ring in rings)
    n_ring_atoms = len({atom for ring in rings for atom in ring})
    if n_ring_atoms == 0:
        return 0.0

    return TRS / n_ring_atoms


class RCICalculator:
    """Callable class returning the Ring Complexity Index."""

    def __call__(self, smiles_or_mol: Union[str, Chem.Mol]) -> float | None:
        return calculate_RCI(smiles_or_mol)

