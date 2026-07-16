[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](https://github.com/AspirinCode/QRCI)
[![ChemRxiv2025](https://img.shields.io/badge/ChemRxiv-10.26434%2Fchemrxiv--2025--mlqwl--v2-green)](https://doi.org/10.26434/chemrxiv-2025-mlqwl-v2)
[![PyPI](https://img.shields.io/badge/PyPI-cyan)](https://pypi.org/project/qrci)

# QRCI

**A Quantitative Ring Complexity Index for Profiling Ring Topology and Chemical Diversity**

![QRCI](https://github.com/AspirinCode/QRCI/blob/main/figures/qrci_cover.png)

## Quantitative Ring Complexity Index

$\mathrm{QRCI}=\frac{\mathrm{TRS}}{N_{\mathrm{ra}}}\left(1+\frac{N_{\mathrm{fr}}}{N_{\mathrm{r}}+1}\right)+\sum_{r}\left[\frac{360}{360-\alpha_{\mathrm{ideal}}(\ell_{r})}\cdot\frac{1}{\ell_{r}}\cdot\lambda_{M}(\ell_{r})\right]+\frac{\sum W_{i}\cdot D_{i}}{\sqrt{N_{\mathrm{ra}}\cdot\mathrm{TRS}}}+\frac{\log(\sqrt{N_{\mathrm{ta}}})}{N_{\mathrm{r}}+1}+W_{m}\cdot\frac{N_{\mathrm{mr}}}{N_{\mathrm{r}}+1}$

* **TRS** (Total Ring Size): Sum of all ring sizes.
* **$N_{\mathrm{ra}}$**: Total number of atoms in all rings.
* **$N_{\mathrm{r}}$**: Total number of rings.
* **$N_{\mathrm{fr}}$** (Fused Rings): Count of rings sharing atoms or bonds.
* **$N_{\mathrm{ta}}$**: Total number of atoms.
* **$N_{\mathrm{mr}}$**: Total number of macrocycles.
* **$W_{m}$**: Weight for macrocycle descriptors.
* **$W_{i}$**: Weight for topological descriptors.
* **$D_{i}$**: Topological ring diversity descriptor.

## Ring Complexity Index

$\mathrm{RCI}(R)=\mathrm{CR}(R)=\frac{\mathrm{SREL}(R)}{\mathrm{SEL}(R)}$

$\mathrm{SREL}(R)$ counts the total number of atom participations across all rings, including duplicate counts if atoms belong to multiple rings. $\mathrm{SEL}(R)$ is the number of unique atoms that appear in at least one ring.

$RCI=\frac{TRS}{nRingAtoms}$, where TRS is the total ring size and `nRingAtoms` is the number of atoms belonging to a ring.

Reference: Gasteiger, J., & Jochum, C. (1979). An Algorithm for the Perception of Synthetically Important Rings. Journal of Chemical Information and Computer Sciences, 19(1), 43-48. https://doi.org/10.1021/ci60017a011

![RCI of Drug](https://github.com/AspirinCode/QRCI/blob/main/figures/drugbank5.1.13_apvd_r1r10w900_rci_histdist.png)

## Installation

The reproducible project environment is defined in [environment.yml](environment.yml):

```bash
conda env create -f environment.yml
conda activate qrci
```

For pip-based local development:

```bash
pip install -r requirements.txt
```

For package-only use from PyPI:

```bash
pip install qrci
```

Current project environment targets:

```text
Python==3.13.2
rdkit==2026.03.4
scipy==1.15.1
```

## Usage

```python
from rdkit import Chem
from QRCI.QRCI import QRCICalculator, get_QRCIproperties
from QRCI.RCI import RCICalculator

smiles = "C1=CCOCc2cc(ccc2OCCN2CCCC2)Nc2nccc(n2)-c2cccc(c2)COC1"

rci_calc = RCICalculator()
print(f"RCI: {rci_calc(smiles):.4f}")

qrci_calc = QRCICalculator(weights="mean")
score_mean = qrci_calc(smiles)
print(f"QRCI(default/mean weights): {score_mean:.4f}")
# QRCI(default/mean weights): 4.0330

mol = Chem.MolFromSmiles(smiles)
props = get_QRCIproperties(mol)
print(props)
```

## Notebooks

* [QRCI_calculate_pip_v2.1.ipynb](QRCI_calculate_pip_v2.1.ipynb): compact QRCI/RCI calculation example using the pip-compatible package API.
* [QRCI_calculate_v1.1.ipynb](QRCI_calculate_v1.1.ipynb): extended analysis workflow with QEPPI, SA score, NP score, QED, and plotting.

## Data

* [DrugBank](https://go.drugbank.com/)
* [iPPI-DB](https://ippidb.pasteur.fr/)
* [COCONUT: the COlleCtion of Open NatUral producTs](https://coconut.naturalproducts.net/)
* [ChEMBL35](https://www.ebi.ac.uk/chembl/)
* [PubChem](https://pubchem.ncbi.nlm.nih.gov/)

## Analysis Helpers

### Molecular Standardization

https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html

https://github.com/rdkit/rdkit/blob/master/Docs/Notebooks/MolStandardize.ipynb

### Spacial Score

```python
from rdkit.Chem.SpacialScore import SPS

sps_score = SPS(mol, normalize=True)
```

https://rdkit.org/docs/source/rdkit.Chem.SpacialScore.html

### SA Score and NP Score

The extended notebook uses RDKit contrib modules:

```python
import os
import sys

sys.path.append(os.path.join(os.environ["CONDA_PREFIX"], "share", "RDKit", "Contrib"))
from SA_Score import sascorer
from NP_Score import npscorer

sascore = sascorer.calculateScore(mol)
fscore = npscorer.readNPModel()
npscore = npscorer.scoreMol(mol, fscore)
```

https://greglandrum.github.io/rdkit-blog/posts/2023-12-01-using_sascore_and_npscore.html

### QED

```python
from rdkit import Chem
from rdkit.Chem import QED

smiles = "C=CCN1CC(C(=O)N(CCCN(C)C)C(=O)NCC)C[C@@H]2c3cccc4[nH]cc(c34)C[C@H]21"
mol = Chem.MolFromSmiles(smiles)

qed_score = QED.qed(mol)
print(f"QED Score: {qed_score:.3f}")
# QED Score: 0.605
```

### QEPPI

QEPPI support is provided by the `QEPPIcommunity` package. The package installs as `QEPPI`:

```python
from rdkit import Chem
from QEPPI import QEPPI_Calculator, get_qeppi_properties

q = QEPPI_Calculator()
q.read()

smiles = "C=CCN1CC(C(=O)N(CCCN(C)C)C(=O)NCC)C[C@@H]2c3cccc4[nH]cc(c34)C[C@H]21"
mol = Chem.MolFromSmiles(smiles)

print(q.qeppi(mol))
print(get_qeppi_properties(mol))
```

[QEPPIcommunity](https://pypi.org/project/QEPPIcommunity/)

https://github.com/AspirinCode/QEPPI-community

## Figures

![Distribution of QRCI for approved drugs](https://github.com/AspirinCode/QRCI/blob/main/figures/drugbank5.1.13_apvd_r1r10w900_qrci_histdist.png)

![Trend of RCI/QRCI over time](https://github.com/AspirinCode/QRCI/blob/main/figures/drug2025_by_year_rci_qrci_Trend_1year_dist.png)

## License

Code is released under the MIT License.

## Cite

* Gasteiger, J. and Jochum, C., 1979. An algorithm for the perception of synthetically important rings. Journal of Chemical Information and Computer Sciences, 19(1), pp.43-48.
* Ertl, P., Schuffenhauer, A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. J Cheminform 1, 8 (2009). https://doi.org/10.1186/1758-2946-1-8
* Krzyzanowski, A., Pahl, A., Grigalunas, M., & Waldmann, H. (2023). Spacial Score: A Comprehensive Topological Indicator for Small-Molecule Complexity. Journal of Medicinal Chemistry, 66(18), 12739-12750. https://doi.org/10.1021/acs.jmedchem.3c00689
* Wang J, Xu K, Ma T, Zhang X, Ma P, Li C, et al. A Quantitative Ring Complexity Index for Profiling Ring Topology and Chemical Diversity. ChemRxiv. 2025; doi: [10.26434/chemrxiv-2025-mlqwl-v2](https://doi.org/10.26434/chemrxiv-2025-mlqwl-v2). This content is a preprint and has not been peer-reviewed.
