# InterCellDB

**InterCellDB** is a tool for analyzing intercellular communication on single-cell RNA sequencing (scRNA-seq) data. It is designed for highly user-defined analysis in specific biological context. To accomplish this,  InterCellDB integrates not only protein-protein interaction database, but also several databases on interaction properties (type, direction, effect) as well as gene product properties (subcellular locations,  molecular functions, involving biological process).

**Key features:**

- Covering human and mouse
- Gene subset selection by subcellular location, molecular function and involving biological process
- Gene pair (inferred from protein-protein interactions) subset selection by type, direction and effect

**Main Functions:**

- Network analysis, which prioritizes significant interactions between some 2-cell groups
- Intercellular analysis on action relations, which shows the composition of action relations between given 2 cell clusters
- Intercellular analysis on filtering gene pairs, which extracts more significant gene pairs for downstream analysis
- Intercellular analysis on specificity, which evaluates co-occurences of gene pair across interactions
- Intercellular analysis on spatial pattern, which shows the summary of target gene pairs with their attributes

## Installation

The installation steps are given as follows: 

# [TODO] check this installation code

``` R
install.packages("devtools")
devtools::install_version(package = "GO.db", version = "3.8.2")
devtools::install_github("ZJUDBlab/InterCellDB")
```

This code is tested to function well in R version 3.6.3. The dependent package `GO.db` is recommended to be version 3.8.2 to accomplish full compatibility.

## Usage of InterCellDB

There are 2 code examples given along with the package. They are located in `code-examples` directory of the package. The example in mouse data shows the routine workflow of InterCellDB, and the example in human data shows the capability of InterCellDB by comparing to other method and the workflow for evaluating off-target effects.

# How to cite

The manuscript is pending.
