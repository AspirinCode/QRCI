from collections import namedtuple
import math
from typing import Iterable, Union

from rdkit import Chem
from rdkit.Chem import Descriptors

from .utils import calculate_strain_factor, count_macrocycles


def calculate_QRCI(
    smiles_or_mol: Union[str, Chem.Mol],
    weights: Iterable[float] | None = None,
    normalize: bool = True,
    W_macro: float = 2.0,
) -> float | None:
    """Compute the Quantitative Ring Complexity Index."""
    if weights is None:
        weights = (1.3, 0.9, 1.2, 0.9, 0.9, 1.1)

    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
        if mol is None:
            return None
    else:
        mol = smiles_or_mol

    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    n_ring = len(rings)
    if n_ring == 0:
        return 0.0

    TRS = sum(len(ring) for ring in rings)
    n_ring_atom = len({atom for ring in rings for atom in ring})
    n_fused_ring = sum(
        1 for ring in rings if any(set(ring) & set(other) for other in rings if ring != other)
    )
    SF = calculate_strain_factor(mol, rings)
    n_total_atom = mol.GetNumAtoms()

    n_arom_hetero = Descriptors.NumAromaticHeterocycles(mol)
    n_arom_carbo = Descriptors.NumAromaticCarbocycles(mol)
    n_ali_hetero = Descriptors.NumAliphaticHeterocycles(mol)
    n_ali_carbo = Descriptors.NumAliphaticCarbocycles(mol)
    n_sat_carbo = Descriptors.NumSaturatedCarbocycles(mol)
    n_sat_hetero = Descriptors.NumSaturatedHeterocycles(mol)
    n_macro = count_macrocycles(mol)

    normalization_factor = (
        math.sqrt(n_ring_atom * TRS) if normalize and n_ring_atom * TRS > 0 else 1
    )

    chemical_diversity = (
        weights[0] * n_arom_hetero
        + weights[1] * n_arom_carbo
        + weights[2] * n_ali_hetero
        + weights[3] * n_ali_carbo
        + weights[4] * n_sat_carbo
        + weights[5] * n_sat_hetero
    ) / normalization_factor

    structural_contribution = (
        TRS / n_ring_atom if n_ring_atom > 0 else 0
    ) * (1 + n_fused_ring / (n_ring + 1))
    molecular_context = (
        math.log(n_total_atom ** 0.5) / (n_ring + 1) if n_total_atom > 1 else 0
    )
    macrocycle_contribution = W_macro * (n_macro / (n_ring + 1)) if n_ring > 0 else 0

    return (
        structural_contribution
        + SF
        + chemical_diversity
        + molecular_context
        + macrocycle_contribution
    )


QRCIproperties = namedtuple(
    "QRCIproperties",
    [
        "nAromHetero",
        "nAromCarbo",
        "nAliHetero",
        "nAliCarbo",
        "nSatHetero",
        "nSatCarbo",
        "nMacrocycles",
        "TRS",
        "nRingAtom",
        "nFusedRing",
        "SF",
    ],
)


def get_QRCIproperties(mol: Chem.Mol) -> QRCIproperties:
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    TRS = sum(len(ring) for ring in rings)
    n_ring_atom = len({atom for ring in rings for atom in ring})
    n_fused_ring = sum(
        1 for ring in rings if any(set(ring) & set(other) for other in rings if ring != other)
    )
    SF = calculate_strain_factor(mol, rings)

    return QRCIproperties(
        nAromHetero=Descriptors.NumAromaticHeterocycles(mol),
        nAromCarbo=Descriptors.NumAromaticCarbocycles(mol),
        nAliHetero=Descriptors.NumAliphaticHeterocycles(mol),
        nAliCarbo=Descriptors.NumAliphaticCarbocycles(mol),
        nSatHetero=Descriptors.NumSaturatedHeterocycles(mol),
        nSatCarbo=Descriptors.NumSaturatedCarbocycles(mol),
        nMacrocycles=count_macrocycles(mol),
        TRS=TRS,
        nRingAtom=n_ring_atom,
        nFusedRing=n_fused_ring,
        SF=SF,
    )


class QRCICalculator:
    def __init__(self, weights: Union[str, Iterable[float]] = "mean", normalize: bool = True, w_macro: float = 2.0):
        self.normalize = normalize
        self.w_macro = w_macro
        self.weights = self._get_weights(weights)

    def _get_weights(self, mode: Union[str, Iterable[float]]):
        if mode == "mean":
            return [1.3, 0.9, 1.2, 0.9, 0.9, 1.1]
        elif mode == "max":
            return [1.5, 1.2, 1.3, 1.1, 1.0, 1.0]
        elif mode == "unit":
            return [1.0] * 6
        elif isinstance(mode, (list, tuple)) and len(mode) == 6:
            return list(mode)
        else:
            raise ValueError(
                "Invalid weight mode or custom weights must be a list/tuple of length 6."
            )

    def __call__(self, mol: Chem.Mol) -> float:
        return calculate_QRCI(
            mol, weights=self.weights, normalize=self.normalize, W_macro=self.w_macro
        )
