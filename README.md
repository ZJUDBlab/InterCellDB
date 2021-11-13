# InterCellDB

**InterCellDB** is a tool for analyzing intercellular communication on single-cell RNA sequencing (scRNA-seq) data. It is designed for highly user-defined analysis in specific biological context. To accomplish this,  InterCellDB integrates not only protein-protein interaction database, but also several databases on interaction properties (type, direction, effect) as well as gene product properties (subcellular locations,  molecular functions, involving biological process).

**Key features:**

- Covering human and mouse
- Gene subset selection by subcellular location, molecular function and involving biological process
- Gene pair (inferred from protein-protein interactions) subset selection by type, direction and effect

**Main Functions:**

- Full network analysis, which summarises the interactions between every 2 cell clusters
- Intercellular analysis on action relations, which shows the composition of action mode and effect between given 2 cell clusters
- Intercellular analysis on filtering gene pairs, which extracts top ranked significant gene pairs for downstream analysis
- Intercellular analysis on specificity, which evaluates co-occurences of gene pair across interactions
- Intercellular analysis on spatial pattern, which shows the summary of target gene pairs with their attributes

## Installation

1. remove the 'GO.db' package 

```R
remove.packages("GO.db")
```

2. install the required packages and re-intall 'GO.db' package from Bioconductor

```R
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install(pkgs = "GO.db")
```

3. install the InterCellDB package from GitHub

```R
devtools::install_github("ZJUDBlab/InterCellDB")
```



## Tutorials

- [Running InterCellDB step-by-step](./tutorials/Basic-steps.md)



# How to cite

The manuscript is pending.
