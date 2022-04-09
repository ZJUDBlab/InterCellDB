# Running InterCellDB step-by-step
## Installation

The installation steps given as follows:

1. remove the 'GO.db' package 

```R
remove.packages("GO.db")
```

2. install the required packages and re-install 'GO.db' package from Bioconductor

```R
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install(pkgs = "GO.db")
```

3. install the InterCellDB package from GitHub

```R
devtools::install_github("ZJUDBlab/InterCellDB")
```



## Running Process


This tutorial illustrates the steps for one routine process using InterCellDB. The main function and features provided by InterCellDB are shown by using the example of the scRNA-seq data from mouse melanoma and lymph node (original article [PMID: 32433953](https://pubmed.ncbi.nlm.nih.gov/32433953/)).

Required libraries for running this tutorial is provided as follows:

```R
# versions of used packages when we tested the code are given
library(InterCellDB)  # 0.9.2
library(Seurat)       # 3.2.0
library(dplyr)        # 1.0.5
library(future)       # 1.21.0
```



## Step 1: Data preparation and create InterCell object

InterCellDB requires 2 data as input:

- normalized count matrix
- differentially expressed genes (DEGs) with their belonging clusters

We provide one pre-processed Seurat object in RDS file to illustrate the data format of the normalized count matrix, and it can be downloaded publicly. 

```R
seurat.obj <- readRDS(url("https://figshare.com/ndownloader/files/31497371"))
tmp.counts <- seurat.obj[["RNA"]]@data
colnames(tmp.counts) <- as.character(Idents(seurat.obj))
# show data structure, row names are genes, column names are clusters
tmp.counts[1:5, 1:5]
# 5 x 5 sparse Matrix of class "dgCMatrix"
#            cDC2        NK        NK       cDC2   Myeloid
# Gnai3 .         1.4526836 0.7893417 0.02189368 .        
# Cdc45 0.1342404 0.1345867 .         1.01913912 0.1698247
# Scml2 .         .         .         .          .        
# Apoh  .         .         .         .          .        
# Narf  .         .         .         0.02189368 .
```

If users are attempting to use their own dataset, make sure the data normalization and cell clustering are done. Many tools can handle with this. Take Seurat for example, the `NormalizeData` or `SCTransform` should be processed for getting normalized data, and `FindClusters` should be processed for getting cell clustering result. 



For DEGs, we also provide the pre-processed RDS file publicly.

```R
markers.used <- readRDS(url("https://figshare.com/ndownloader/files/31497401"))
#                p_val    LogFC pct.1 pct.2          PVal cluster    gene
# Plbd1   2.436468e-234 1.751750 0.990 0.293 4.703844e-230    cDC2   Plbd1
# Cd209a  1.126670e-219 2.934390 0.642 0.099 2.175148e-215    cDC2  Cd209a
# Ms4a6c  2.003678e-217 1.306534 0.916 0.244 3.868301e-213    cDC2  Ms4a6c
# Alox5ap 2.043584e-213 1.118439 0.922 0.214 3.945343e-209    cDC2 Alox5ap
# Fgr     3.582238e-205 1.177831 0.775 0.154 6.915869e-201    cDC2     Fgr
```

To note, InterCellDB defines several mandatory column names when read in DEGs, which are columns named 'gene', 'cluster', 'LogFC', 'PVal', and their meanings are listed below:

- Column 'gene', the gene name.
- Column 'cluster', the belonging cluster of every DEG.
- Column 'LogFC', the log fold changes of each DEG.
- Column 'PVal', the statistical significance for the gene being defined as DEG.



Then, the `InterCell` object is created and we name it here as `inter.obj`. This variable is used for doing all the following analysis.
```R
inter.obj <- CreateInterCellObject(markers.used, species = "mouse", add.exprs = TRUE, exprs.data = tmp.counts)
```




## Step 2: Customize interactions and proteins from databases

InterCellDB provides customized settings in selecting interactions and proteins.

For customizing interactions, InterCellDB provides 4 options: evidence source, credibility score, action mode and action effect. To run this example data, we use all experimentally validated interactions, pathway-curated interactions and predicted interactions, and select those physically associated ones.

```R
inter.obj <- SelectDBSubset(inter.obj,
		use.exp = TRUE,   # to use experimentally validated interactions
		exp.score.range = c(1, 1000),  # use all credibility score for 'use.exp'
		use.know = TRUE,  # to use pathway curated interactions
		know.score.range = c(1, 1000),  # use all credibility score for 'use.know'
		use.pred = TRUE,  # to use predicted interactions
		pred.score.range = c(1, 1000),  # use all credibility score for 'use.pred' 
		sel.action.mode = "binding",  # select physically associated ones
		sel.action.effect = "ALL")    # not limiting action effects
```


Details about 4 options on customizing interaction database are given in the supporting information of the article. Here, we list some of the most commonly used options.

- evidence sources: experimentally validated (`use.exp`), pathway curated (`use.know`), predicted (`use.pred`).
- credibility score: ranging from 1 to 1000. Highest confidence should >= 900, and high confidence should >= 700, and medium >= 400, while low < 400.
- action mode: giving the relation of interacting proteins. Commonly used one is 'binding', which selects physically associated protein interactions.
- action effect: giving the effect of protein interaction. Commonly used ones are 'positive' and 'negative'. 'positive' means one protein promotes the expression of the other, while 'negative' means the inhibition of the other.



For customizing included proteins, InterCellDB also provides 4 options: protein expression change, protein subcellular location, function and relevant biological process. Here, we select 2 subsets of proteins using InterCellDB on the latter 3 options and fetch their corresponding gene symbols. One subset is the genes for signal receiving cells (denoted as `genes.receiver` below) and the other for signal sending cells (denoted as `genes.sender` below).

```R
# fetch genes of interest. GO.db version is v3.8.2 here, different version may have non-identical result.
genes.receiver <- FetchGeneOI(inter.obj, 
		sel.location = "Plasma Membrane",  # fetch proteins located in plasma membrane
		sel.location.score = c(4:5),       # 4 and 5 are of high confidence
		sel.type = "Receptor",             # fetch those receptors
		sel.go.terms = "GO:0006955"        # fetch [GO:0006955]immune response related proteins
		)
# [1] "Fetch 160 genes of interest."
genes.sender <- FetchGeneOI(inter.obj, 
		sel.location = "Extracellular Region",  # fetch proteins located in extracellular region
		sel.location.score = c(4:5),            # 4 and 5 are of high confidence
		sel.go.terms = "GO:0006955"             # fetch [GO:0006955]immune response related proteins
		)
# [1] "Fetch 281 genes of interest."
```

To note, as the protein expression change is highly associated with the input data, InterCellDB put it to be done in later steps, see parameter `sel.exprs.change` in  `AnalyzeInterInFullView`. If users want to explore more selections on the other 3 options, InterCellDB provides several functions to list all items for them.

- protein subcellular location: `ListAllGeneLocation`
- protein function: `ListAllGeneType` and `ListAllGeneMergeType`
- protein involved biological process: `ListAllGeneGOTerm`



## Step 3: Perform full network analysis

Full network analysis summarizes the interactions between every 2 cell clusters, and infers main participants by aggregated power and total count of gene pairs. The power of every gene pair is calculated by the product of expressions of 2 participating genes. The selected gene pairs are those with *p*-value < 0.05 in cell label permutation test (the *p*-value cutoff can also be set by users).

Full network analysis goes with 2 steps:

1. generate cell label permutation list
2. network analysis and filter statistically significant gene pairs

Here, we do cell label permutation on normalized count matrix for 100 times with 2 parallel processes. During every time of permutation, the average gene expressions for each gene and for every cell cluster are re-calculated. 

```R
plan("multiprocess", workers = 2)  # package future provides the parallel interface, here create 2 parallel processes
tmp.permlist <- Tool.GenPermutation(inter.obj, tmp.counts, perm.times = 100)
plan("sequential")  # close the parallel processes
```

Then, we perform the full network analysis using the previously customized genes and interactions with 4 parallel processes. The one-tailed statistical test is performed with *p*-value < 0.05. 

```R
plan("multiprocess", workers = 4)
inter.obj <- AnalyzeInterInFullView(inter.obj, 
		sel.some.genes.X = genes.receiver,  # set the genes for receiver cells
		sel.some.genes.Y = genes.sender,    # set the genes for sender cells
		sel.exprs.change = "Xup.Yup",       # select genes up-regulated for both sender cells and receiver cells
		run.permutation = TRUE,             # run statistical test with permutation list
		perm.expression = tmp.permlist,     # given the permutation list, generated by `Tool.GenPermutation`
		perm.pval.sides = "one-side-large", # use one-tailed test
		perm.pval.cutoff = 0.05)            # p-value cutoff 0.05
plan("sequential")
```

The results of full network analysis can be further collected by `GetResultFullView` for visualization or table output. Users can define their clusters of interest. Here, we select the interactions among clusters in tumor from example data, which exclude interactions involving 'LN T cell', 'Endo lymphatic', 'Endo LN' and 'Fibroblast LN' clusters.

```R
torm.LN.clusters <- c("LN T cell", "Endo lymphatic", "Endo LN", "Fibroblast LN")
all.clusters <- ListAllClusters(inter.obj)                # list all clusters
used.clusters <- setdiff(all.clusters, torm.LN.clusters)  # select those clusters in tumor
used.clusters <- used.clusters[order(used.clusters)]
# show the result of full network analysis
fullview.result <- GetResultFullView(inter.obj, 
		show.clusters.in.x = used.clusters,
		show.clusters.in.y = used.clusters,
		plot.axis.x.name = "signal receiving cells",
		plot.axis.y.name = "signal sending cells",
		nodes.size.range = c(1,8))
Tool.ShowGraph(fullview.result)  # visualization
Tool.WriteTables(fullview.result, dir.path = "./")  # write result into csv file
```

<img src="./files-Basic-steps/fullview-res.png" alt="fullview-res" style="zoom:25%; align-items:center;" />



## Step 4: Perform intercellular analysis

Intercellular analysis focuses on one cell-cell interaction, and explores the results in 4 aspects:

1. summarize action mode and action effect
2. subset and rank gene pairs
3. evaluate the specificity of gene pairs
4. summarize gene pairs in spatial pattern



With the help of full network analysis, we can pick up one highly active interaction or collect one interaction of our interest. Here, we choose the interaction between Myeloid and CAF 1. CAF 1 is treated as the signal sending cell (in Y axis of the graph above), and Myeloid is treated as the signal receiving cells (in X axis). 

```R
inter.obj <- FetchInterOI(inter.obj, cluster.x = "Myeloid", cluster.y = "CAF 1")
```



### 4.1 Summarize action mode and action effect

The summary on action mode and action effect can be generated using `AnalyzeInterInAction`. 

```R
inter.obj <- AnalyzeInterInAction(inter.obj)
```

Then, we can show the composition of action mode for all gene pairs. Here, we use 6 action modes.

```R
used.action.mode <- c("activation", "inhibition", "catalysis", "reaction", "expression", "ptmod")  # action mode: not showing mode 'binding' and 'other'
used.color.mode <- c("#FB8072", "#80B1D3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDB462")  # customized color, one-to-one corresponding to those in `used.action.mode`
action.mode.result <- GetResultPieActionMode(inter.obj, 
		limits.exprs.change = c("Xup.Yup"),
		limits.action.mode = used.action.mode,
		color.action.mode = used.color.mode)
Tool.ShowGraph(action.mode.result)
```
<img src="./files-Basic-steps/act-mode.png" alt="act-mode" style="zoom:20%; align-items:center;" />

We can also show the composition of action effect.

```R
action.effect.result <- GetResultPieActionEffect(inter.obj, limits.exprs.change = c("Xup.Yup"))
Tool.ShowGraph(action.effect.result)
```
<img src="./files-Basic-steps/act-effect.png" alt="act-effect" style="zoom:20%; align-items:center;" />



### 4.2 Subset and rank gene pairs

Circumstances exist when too many gene pairs are returned or the initial selection of gene subsets are not well satisfied to users' needs. We provide additional function `SelectInterSubset` for picking up subsets of gene pairs in intercellular scale. All subset selection given in full network analysis are applicable in this function. In addition, action mode and action effect can be customized. 

Here, we filter activated gene pairs, i.e. selecting those with either 'activation' action mode or 'positive' action effect.

```R
inter.obj <- SelectInterSubset(inter.obj, 
		sel.action.mode = "activation",
		sel.action.effect = "positive",
		sel.action.merge.option = "union"
		)
```

Then, we visualize all gene pairs with their power and confidence. 

```R
result.inter.pairs <- GetResultTgCrosstalk(inter.obj, 
		func.to.use = list(
			Power = function(x) { log1p(x) },        # further log-transform the power
			Confidence = function(x) { 1/(1+x) }     # invert the confidence to show
		),  
		axis.order.xy = c("Power", "Power"),       # order genes by power
		axis.order.xy.decreasing = c(TRUE, TRUE),  # select if ordering is decreasing
		nodes.size.range = c(1, 8))
Tool.ShowGraph(result.inter.pairs)
```

<img src="./files-Basic-steps/inter-rank-gene-pair.png" alt="inter-rank-gene-pair" style="zoom:20%; align-items:center;" />

We rank gene pairs by their power, and select those top ranked genes.

```R
inter.obj <- SelectInterSubset(inter.obj, 
		sel.some.genes.X = c("Itgb2","Itgam","C5ar1","Ccr2","C3ar1","Ccr1","Ccr5"),
		sel.some.genes.Y = c("C3","Ccl11","Cxcl12","Cxcl1","Ccl2","Ccl7","C4b")
		)
```



### 4.3 Evaluate the specificity of gene pairs

Besides the power rank, we also evaluate the specificity of the gene pairs, i.e. checking whether one gene pair is exclusively presented in one interaction or widely observed. 

Here, we choose all clusters in tumor as the sender cells and 'Myeloid' cluster as the receiver cell. This is to evaluate specificity of the signaling protein that 'Myeloid' cell receives.

```R
tmp.target.cluster.groups <- ListClusterGroups(inter.obj, 
		use.former = TRUE,
		cluster.former = c("Myeloid")  # list all crosstalk with Myeloid as receiver cell, this order is aligned with the X and Y axis in fullview result. In this example, clusters listed in X axis are the receiver cells.
		)
tmp.target.cluster.groups <- setdiff(tmp.target.cluster.groups, c("Myeloid~Myeloid", "Myeloid~Endo LN", "Myeloid~Endo lymphatic", "Myeloid~Fibroblast LN", "Myeloid~LN T cell"))
```

The analysis is processed in `AnalyzeInterSpecificity` and users can visualize the result by `GetResultTgSpecificity`.

```R
inter.obj <- AnalyzeInterSpecificity(inter.obj, 
		to.cmp.cluster.groups = tmp.target.cluster.groups
		)

result.specificity <- GetResultTgSpecificity(inter.obj,
		sel.uq.cnt.options = seq_along(tmp.target.cluster.groups),
		plot.uq.cnt.merged = TRUE, 
		plot.name.X.to.Y = FALSE,
		func.to.use = list(Power = function(x) { log1p(x) },
                       Confidence = function(x) { 1/(1+x) }),
		dot.size.range = c(2,8),
		dot.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
		dot.colour.value.seq = c(0, 0.4, 1)
		)
Tool.ShowGraph(result.specificity)
```

<img src="./files-Basic-steps/eval-specificity.png" alt="eval-specificity" style="zoom:20%; align-items:center;" />



### 4.4 Summarize gene pairs in spatial pattern

The protein-protein interaction is limited by their spatial closeness. InterCellDB provides genes with their corresponding gene product subcellular locations, and also provides visualization method to show gene pairs in spatial pattern. The action mode, action effect, protein subcellular location, and gene expression change can be integrated into the plot. 

Here, we choose only the plasma membrane protein and extracellular region protein to show the paracrine protein-protein interactions.

```R
tmp.hide.locations <- setdiff(ListAllGeneLocation(inter.obj), c("Plasma Membrane", "Extracellular Region"))  # select the subcellular locations allowed to show
set.seed(101L)  # set seed to make reproducible result, or gene will be re-arranged every time re-running GetResultTgCellPlot
result.sptialpattern <- GetResultTgCellPlot(inter.obj, 
		area.extend.times = 20,  # control the size of plotting. Increase the value by 10 when warnings ask to
		hide.other.area = TRUE,
		hide.locations.X = tmp.hide.locations,
		hide.locations.Y = tmp.hide.locations,
    expand.gap.radius.list = list(ECM = 8, CTP = 2, NC = 2, OTHER = 2),
		link.size = 0.3,
		link.alpha = 0.8,
		legend.manual.left.spacing = grid::unit(0.1, "cm")
		)
Tool.ShowGraph(result.sptialpattern)
```

<img src="./files-Basic-steps/spatial-pattern.png" alt="spatial-pattern-plot" style="zoom:20%; align-items:center;" />







**Session Info**

```R
sessionInfo()
```

```R
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] future_1.21.0     dplyr_1.0.5       Seurat_3.2.0      InterCellDB_0.9.2

loaded via a namespace (and not attached):
  [1] Rtsne_0.15           colorspace_2.0-0     deldir_0.2-10       
  [4] ellipsis_0.3.1       ggridges_0.5.3       spatstat.data_2.1-0 
  [7] farver_2.1.0         leiden_0.3.7         listenv_0.8.0       
 [10] ggrepel_0.9.1        bit64_4.0.5          AnnotationDbi_1.46.1
 [13] fansi_0.4.2          codetools_0.2-18     splines_3.6.3       
 [16] cachem_1.0.4         polyclip_1.10-0      jsonlite_1.7.2      
 [19] ica_1.0-2            cluster_2.1.1        GO.db_3.8.2         
 [22] png_0.1-7            uwot_0.1.10          shiny_1.6.0         
 [25] sctransform_0.3.2    compiler_3.6.3       httr_1.4.2          
 [28] assertthat_0.2.1     Matrix_1.3-2         fastmap_1.1.0       
 [31] lazyeval_0.2.2       later_1.1.0.1        htmltools_0.5.1.1   
 [34] prettyunits_1.1.1    tools_3.6.3          rsvd_1.0.3          
 [37] igraph_1.2.6         gtable_0.3.0         glue_1.4.2          
 [40] RANN_2.6.1           reshape2_1.4.4       Rcpp_1.0.7          
 [43] spatstat_1.64-1      Biobase_2.44.0       vctrs_0.3.7         
 [46] ape_5.4-1            nlme_3.1-152         lmtest_0.9-38       
 [49] stringr_1.4.0        globals_0.14.0       mime_0.10           
 [52] miniUI_0.1.1.1       lifecycle_1.0.0      irlba_2.3.3         
 [55] goftest_1.2-2        MASS_7.3-53.1        zoo_1.8-9           
 [58] scales_1.1.1         hms_1.0.0            promises_1.2.0.1    
 [61] spatstat.utils_2.1-0 parallel_3.6.3       RColorBrewer_1.1-2  
 [64] memoise_2.0.0        reticulate_1.18      pbapply_1.4-3       
 [67] gridExtra_2.3        ggplot2_3.3.3        rpart_4.1-15        
 [70] stringi_1.5.3        RSQLite_2.2.5        S4Vectors_0.22.1    
 [73] BiocGenerics_0.30.0  rlang_0.4.10         pkgconfig_2.0.3     
 [76] matrixStats_0.58.0   lattice_0.20-41      ROCR_1.0-11         
 [79] purrr_0.3.4          tensor_1.5           labeling_0.4.2      
 [82] patchwork_1.1.1      htmlwidgets_1.5.3    cowplot_1.1.1       
 [85] bit_4.0.4            tidyselect_1.1.0     parallelly_1.24.0   
 [88] RcppAnnoy_0.0.18     plyr_1.8.6           magrittr_2.0.1      
 [91] R6_2.5.0             IRanges_2.18.3       generics_0.1.0      
 [94] DBI_1.1.1            mgcv_1.8-34          pillar_1.5.1        
 [97] fitdistrplus_1.1-3   survival_3.2-10      abind_1.4-5         
[100] tibble_3.1.0         future.apply_1.7.0   crayon_1.4.1        
[103] KernSmooth_2.23-18   utf8_1.2.1           plotly_4.9.3        
[106] progress_1.2.2       grid_3.6.3           data.table_1.13.6   
[109] blob_1.2.1           digest_0.6.27        xtable_1.8-4        
[112] tidyr_1.1.3          httpuv_1.5.5         stats4_3.6.3        
[115] munsell_0.5.0        viridisLite_0.3.0   
```

