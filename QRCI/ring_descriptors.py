"""Ring descriptor helper functions for QRCI package."""

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, rdMolDescriptors


# Calculates the total number of stereocenters
def get_num_stereocenters(mol: Chem.Mol) -> int:
    """Return the total number of specified and unspecified stereocenters."""
    return AllChem.CalcNumAtomStereoCenters(mol) + AllChem.CalcNumUnspecifiedAtomStereoCenters(mol)


# calculate the number of macrocycle
def get_num_macrocycles(mol: Chem.Mol, min_ring_size: int = 12) -> int:
    """Return the number of rings with size >= ``min_ring_size``."""
    if mol is None:
        return 0
    sssr = Chem.GetSymmSSSR(mol)
    count = 0
    for ring in sssr:
        if len(ring) >= min_ring_size:
            count += 1
    return count


# Calculate the number of heteroatoms
def calculate_number_of_heteroatoms(mol: Chem.Mol) -> int | None:
    """Return the number of heteroatoms (non H or C) in the molecule."""
    if mol is None or not hasattr(mol, "GetAtoms"):
        return None
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    return len(heteroatoms)


# calculate the number of ring atoms
def count_ring_atoms(mol: Chem.Mol) -> int | None:
    """Return the number of atoms that are part of any ring."""
    if mol is None:
        return None
    return sum(1 for atom in mol.GetAtoms() if atom.IsInRing())


# Calculates the ratio of heteroatoms to heavy atoms
def calculate_heteroatom_ratio(mol: Chem.Mol) -> float | None:
    """Return the ratio of heteroatoms to heavy atoms."""
    if mol is None:
        return None
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms == 0:
        return 0.0
    num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    heteroatom_ratio = num_heteroatoms / num_heavy_atoms
    return heteroatom_ratio


# calculate the number of isolated rings
def calculate_nIR(mol: Chem.Mol) -> int | None:
    """Return the number of isolated rings (nIR)."""
    if mol is None:
        return None
    sssr = rdmolops.GetSymmSSSR(mol)
    ring_list = [set(ring) for ring in sssr]
    isolated_rings = 0
    for i, ring_i in enumerate(ring_list):
        is_isolated = all(ring_i.isdisjoint(ring_j) for j, ring_j in enumerate(ring_list) if i != j)
        if is_isolated:
            isolated_rings += 1
    return isolated_rings


# calculate the number of fused rings
def calculate_nFR(mol: Chem.Mol) -> int | None:
    """Return the number of fused rings (nFR)."""
    if mol is None:
        return None
    sssr = rdmolops.GetSymmSSSR(mol)
    total_rings = len(sssr)
    nIR = calculate_nIR(mol)
    nFR = total_rings - nIR
    return nFR


# calculate the total ring size (TRS)
def calculate_TRS(mol: Chem.Mol) -> int | None:
    """Return the total ring size (sum of all ring sizes)."""
    if mol is None:
        return None
    sssr = rdmolops.GetSymmSSSR(mol)
    total_ring_size = sum(len(ring) for ring in sssr)
    return total_ring_size


# calculate sum of all bonds forming the rings
def calculate_Rperim(mol: Chem.Mol) -> int | None:
    """Return the total ring perimeter (equivalent to ``calculate_TRS``)."""
    if mol is None:
        return None
    return calculate_TRS(mol)


# Calculate the Number of Ring Systems (NRS)
def calculate_nrs(mol: Chem.Mol) -> int | None:
    """Return the number of ring systems (NRS)."""
    if mol is None:
        return None
    nBO = mol.GetNumBonds()
    nSK = mol.GetNumHeavyAtoms()
    B_R = sum(bond.IsInRing() for bond in mol.GetBonds())
    A_R = sum(atom.IsInRing() for atom in mol.GetAtoms())
    nStructures = len(rdmolops.GetMolFrags(mol, asMols=True))
    NRS = (nBO - B_R) - (nSK - A_R) + nStructures
    return NRS


# calculate Normalized Number of Ring Systems (NNRS)
def calculate_cyclomatic_number(mol: Chem.Mol) -> int | None:
    """Return the cyclomatic number (nCIC)."""
    if mol is None:
        return None
    nBO = mol.GetNumBonds()
    nSK = mol.GetNumHeavyAtoms()
    nStructures = len(rdmolops.GetMolFrags(mol, asMols=True))
    nCIC = nBO - nSK + nStructures
    return nCIC


def calculate_nnrs(mol: Chem.Mol) -> float | None:
    """Return the Normalized Number of Ring Systems (NNRS)."""
    if mol is None:
        return None
    nrs = calculate_nrs(mol)
    nCIC = calculate_cyclomatic_number(mol)
    if nCIC is None or nCIC == 0:
        return None
    nnrs = nrs / nCIC
    return nnrs


# Molecular Cyclized Degree (MCD)
def calculate_mcd(mol: Chem.Mol) -> float | None:
    """Return the Molecular Cyclized Degree (MCD)."""
    if mol is None:
        return None
    nSK = mol.GetNumHeavyAtoms()
    nRingAtoms = sum(atom.IsInRing() for atom in mol.GetAtoms())
    if nSK == 0:
        return None
    MCD = nRingAtoms / nSK
    return MCD


# calculate Ring Fusion Density (RFD)
def calculate_ring_fusion_density(mol: Chem.Mol) -> float | None:
    """Return the Ring Fusion Density (RFD)."""
    if mol is None:
        return None
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    total_rings = len(atom_rings)
    if total_rings == 0:
        return None
    fused_rings = 0
    for i, ring1 in enumerate(atom_rings):
        for j, ring2 in enumerate(atom_rings):
            if i < j and set(ring1).intersection(set(ring2)):
                fused_rings += 1
                break
    rfd = fused_rings / total_rings
    return rfd


# calculate Aromatic Ratio (ARR)
def Calc_ARR(mol: Chem.Mol) -> float:
    """Return the Aromatic Ratio (ARR)."""
    m = Chem.RemoveHs(mol)
    num_bonds = m.GetNumBonds()
    if num_bonds == 0:
        return 0.0
    num_aromatic_bonds = sum(1 for bond in m.GetBonds() if bond.GetIsAromatic())
    ARR = num_aromatic_bonds / num_bonds
    return round(ARR, 4)


# Calculate the aromatic-sp3 carbon balance for a molecule.
def Calc_Ar_Alk_balance(mol: Chem.Mol) -> int:
    """Return difference between aromatic and sp3 carbons."""
    m = Chem.RemoveHs(mol)
    num_aromatic_carbon = sum(
        1 for atom in m.GetAromaticAtoms() if atom.GetSymbol() == "C"
    )
    num_sp3_carbon = sum(
        1
        for atom in m.GetAtoms()
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and atom.GetSymbol() == "C"
    )
    return num_aromatic_carbon - num_sp3_carbon


# Calculates the ratio of HeteroRings to all rings
def calculate_heterorings_ratio(mol: Chem.Mol) -> float | None:
    """Return the ratio of heteroatom-containing rings to total rings."""
    if mol is None:
        return None
    total_rings = rdMolDescriptors.CalcNumRings(mol)
    if total_rings == 0:
        return 0.0
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    hetero_rings = sum(
        any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in (6, 1) for idx in ring)
        for ring in atom_rings
    )
    heterorings_ratio = hetero_rings / total_rings
    return round(heterorings_ratio, 4)
