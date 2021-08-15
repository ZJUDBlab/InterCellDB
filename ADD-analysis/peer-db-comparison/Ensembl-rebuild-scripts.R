# This is the scripts for ensembl database rebuilt.

entrez.db <- entrez.9606.db


# get new gene stable ID
new.gene.id <- "ENSG00000205777"
new.gene.id
ensembl.g37.db[which(ensembl.g37.db$Gene.stable.ID == new.gene.id), ]
ensembl.g38.db[which(ensembl.g38.db$Gene.stable.ID == new.gene.id), ]

new.gene.name <- "GAGE1"
new.gene.name
entrez.db[which(entrez.db$Symbol_from_nomenclature_authority == new.gene.name), ]

# replace
t.i <- 3
test.check[t.i, ]
#
test.check[t.i, "new.Gene.stable.ID"] <- 
test.check[t.i, "new.Gene.name"] <- 
#
test.check[(t.i-2):(t.i+1),]