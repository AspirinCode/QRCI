# QRCI

**Quantitative Ring Complexity Index for Evaluating Molecular Structure and Chemical Diversity** 

![QRCI](https://github.com/AspirinCode/qrci/blob/main/figures/qrci_cover.png)

### Quantitative Ring Complexity Index

$\mathrm{QRCI}=\frac{\mathrm{TRS}}{N_{\mathrm{ra}}}\left(1+\frac{N_{\mathrm{fr}}}{N_{\mathrm{r}}+1}\right)+\sum_{r}\left[\frac{360}{360-\alpha_{\mathrm{ideal}}(\ell_{r})}\cdot\frac{1}{\ell_{r}}\cdot\lambda_{M}(\ell_{r})\right]+\frac{\sum W_{i}\cdot D_{i}}{\sqrt{N_{\mathrm{ra}}\cdot\mathrm{TRS}}}+\frac{\log(N_{\mathrm{ta}})}{N_{\mathrm{r}}+1}+W_{m}\cdot\frac{N_{\mathrm{mr}}}{N_{\mathrm{r}}+1}$  

* **TRS** (Total Ring Size): Sum of all ring sizes.
* **$N_{\mathrm{ra}}$**: Total number of atoms in all rings.
* **$N_{\mathrm{r}}$**: Total number of rings
* **$N_{\mathrm{fr}}$** (Fused Rings): Count of rings sharing atoms or bonds.
* **$N_{\mathrm{ta}}$**: Total number of atoms
* **$N_{\mathrm{mr}}$**: total number of macrocycles
* **$W_{m}$**: Weight for macrocycle descriptors.
* **$W_{i}$**: Weight for topological descriptors.
* **$D_{i}$**: Topological ring diversity descriptor.


### Ring Complexity Index
$RCI=\frac{TRS}{nRingAtoms}$  
where TRS is the total ring size and $nRingAtoms$ is the number of atoms belonging to a ring.
Ref: Gasteiger, J., & Jochum, C. (1979). An Algorithm for the Perception of Synthetically Important Rings. Journal of Chemical Information and Computer Sciences, 19(1), 43–48. https://doi.org/10.1021/ci60017a011  


## Requirements
```python
Python==3.13.2
rdkit==2024.09.6
scipy==1.15.1
```

## Data

* [DrugBank](https://go.drugbank.com/)  
* [iPPI-DB](https://ippidb.pasteur.fr/)  
* [COCONUT: the COlleCtion of Open NatUral producTs](https://coconut.naturalproducts.net/)  
* [ChEMBL35](https://www.ebi.ac.uk/chembl/)  


### Molecular Standardization

https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html

https://github.com/rdkit/rdkit/blob/master/Docs/Notebooks/MolStandardize.ipynb



## QRCI calculation

```python
# calculate QRCI with normalization and entropy-based weight optimization
calculate_QRCI(smiles_or_mol, weights=None, normalize=True, W_macro=2.0)
```


## Analysis

### Spacial Score

```python
rdkit.Chem.SpacialScore.SPS(mol, normalize=True)
```

https://rdkit.org/docs/source/rdkit.Chem.SpacialScore.html  

https://github.com/frog2000/Spacial-Score  


### SAscore

```python
#Calculating NP_Score
import npscorer
fscore = npscorer.readNPModel()
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
#QED Score: 0.605
```

https://www.rdkit.org/docs/source/rdkit.Chem.QED.html#module-rdkit.Chem.QED


### QEPPI
quantitative estimate of protein-protein interaction targeting drug-likeness  

```python
#Calculates QEPPI
q = ppi.QEPPI_Calculator()
print("QEPPI.model LOADING...")
q.load("./QEPPI/QEPPI.model")
```

https://github.com/ohuelab/QEPPI  



### Others

**Drug Data From the ChEMBL**

https://github.com/PatWalters/practical_cheminformatics_tutorials/tree/main/misc


## License
Code is released under MIT LICENSE.


## Cite

* Gasteiger, J. and Jochum, C., 1979. An algorithm for the perception of synthetically important rings. Journal of Chemical Information and Computer Sciences, 19(1), pp.43-48.
* Ertl, P., Schuffenhauer, A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. J Cheminform 1, 8 (2009). https://doi.org/10.1186/1758-2946-1-8
* Krzyzanowski, A., Pahl, A., Grigalunas, M., & Waldmann, H. (2023). Spacial Score─A Comprehensive Topological Indicator for Small-Molecule Complexity. Journal of medicinal chemistry, 66(18), 12739–12750. https://doi.org/10.1021/acs.jmedchem.3c00689  





