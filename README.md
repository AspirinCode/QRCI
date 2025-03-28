# QRCI

**Quantitative Ring Complexity Index for Evaluating Molecular Structure and Chemical Diversity** 

![QRCI](https://github.com/AspirinCode/qrci/blob/main/figures/qrci_cover.png)

### Quantitative Ring Complexity Index

$\mathrm{QRCI}=\frac{\mathrm{TRS}}{N_{\mathrm{ra}}}\left(1+\frac{N_{\mathrm{fr}}}{N_{\mathrm{r}}+1}\right)+\sum_{r}\left[\frac{360}{360-\alpha_{\mathrm{ideal}}(\ell_{r})}\cdot\frac{1}{\ell_{r}}\cdot\lambda_{M}(\ell_{r})\right]+\frac{\sum W_{i}\cdot D_{i}}{\sqrt{N_{\mathrm{ra}}\cdot\mathrm{TRS}}}+\frac{\log(N_{\mathrm{ta}})}{N_{\mathrm{r}}+1}+W_{m}\cdot\frac{N_{\mathrm{mr}}}{N_{\mathrm{r}}+1}$  

* **TRS** (Total Ring Size): Sum of all ring sizes.
* **nRingAtom**: Total number of atoms in all rings.
* **nRing**: Total number of rings
* **nFR** (Fused Rings): Count of rings sharing atoms or bonds.
* **nTotalAtom**: Total number of atoms
* **nMacrocycle**: total number of macrocycles



### Ring Complexity Index
$RCI=\frac{TRS}{nRingAtoms}$  
where TRS is the total ring size and $nRingAtoms$ is the number of atoms belonging to a ring.
Ref: Gasteiger, J., & Jochum, C. (1979). An Algorithm for the Perception of Synthetically Important Rings. Journal of Chemical Information and Computer Sciences, 19(1), 43–48. https://doi.org/10.1021/ci60017a011  


## Requirements
```python
Python==3.13.1
rdkit==2024.09.6
scipy==1.14.1
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


```


## Analysis

### Spacial Score

```python
rdkit.Chem.SpacialScore.SPS(mol, normalize=True)
```

https://rdkit.org/docs/source/rdkit.Chem.SpacialScore.html  

https://github.com/frog2000/Spacial-Score  



## License
Code is released under MIT LICENSE.


## Cite

* Gasteiger, J. and Jochum, C., 1979. An algorithm for the perception of synthetically important rings. Journal of Chemical Information and Computer Sciences, 19(1), pp.43-48.
* Krzyzanowski, A., Pahl, A., Grigalunas, M., & Waldmann, H. (2023). Spacial Score─A Comprehensive Topological Indicator for Small-Molecule Complexity. Journal of medicinal chemistry, 66(18), 12739–12750. https://doi.org/10.1021/acs.jmedchem.3c00689  





