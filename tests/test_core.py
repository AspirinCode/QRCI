from rdkit import Chem

from QRCI.QRCI import QRCICalculator, calculate_QRCI, get_QRCIproperties
from QRCI.RCI import RCICalculator, calculate_RCI


PACRITINIB = "C1=CCOCc2cc(ccc2OCCN2CCCC2)Nc2nccc(n2)-c2cccc(c2)COC1"


def test_rci_accepts_smiles_and_mol():
    mol = Chem.MolFromSmiles(PACRITINIB)

    assert calculate_RCI(PACRITINIB) == calculate_RCI(mol)
    assert RCICalculator()(PACRITINIB) > 0


def test_qrci_accepts_smiles_and_mol():
    mol = Chem.MolFromSmiles(PACRITINIB)
    calc = QRCICalculator(weights="mean")

    assert calculate_QRCI(PACRITINIB) == calculate_QRCI(mol)
    assert calc(PACRITINIB) > 0


def test_qrci_properties_include_ring_counts():
    mol = Chem.MolFromSmiles(PACRITINIB)
    props = get_QRCIproperties(mol)

    assert props.TRS > 0
    assert props.nRingAtom > 0
    assert props.nMacrocycles >= 0


def test_invalid_smiles_returns_none():
    assert calculate_RCI("not-a-smiles") is None
    assert calculate_QRCI("not-a-smiles") is None