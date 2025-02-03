# QRCI
**Quantitative Ring Complexity Index for Evaluating Molecular Structure and Chemical Diversity** 

### Quantitative Ring Complexity Index

QRCI=\frac{TRS}{nRingAtoms}\cdot\left(1+\frac{nFusedRings}{nRings+1}\right)+SF+\frac{\sum W_i\cdot D_i}{\sqrt{nRingAtoms\cdot TRS}}+\frac{\log(nTotalAtoms)}{nRings+1}+W_{macro}\cdot\frac{nMacrocycles}{nRings+1}


### Ring Complexity Index
$RCI=\frac{TRS}{nRingAtoms}$  
where TRS is the total ring size and $nRingAtoms$ is the number of atoms belonging to a ring.
Ref: Gasteiger, J., & Jochum, C. (1979). An Algorithm for the Perception of Synthetically Important Rings. Journal of Chemical Information and Computer Sciences, 19(1), 43–48. https://doi.org/10.1021/ci60017a011  


## Requirements
```python
Python==3.13
rdkit==2024.09.5

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



## License
Code is released under MIT LICENSE.


## Cite

* Gasteiger, J. and Jochum, C., 1979. An algorithm for the perception of synthetically important rings. Journal of Chemical Information and Computer Sciences, 19(1), pp.43-48.





