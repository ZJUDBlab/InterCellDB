
# This file is to check data given in Figure 1


#### Human part
# ensembl.9606.db get 127750 proteins (inferred from protein.stable.ID)

human.orig.protein <- unique(c(string.9606.readin.db[,1], string.9606.readin.db[,2]))

human.split.protein <- lapply(human.orig.protein, function(x) {strsplit(x, split=".", fixed=TRUE) })
human.split.protein <- unlist(human.split.protein)  
human.split.protein <- human.split.protein[2* (1:(length(human.split.protein)/2))]
human.split.protein <- unique(human.split.protein)


to.use.human.orig <- unique(intersect(human.split.protein, ensembl.9606.db[,2]))
length(to.use.human.orig) # 19327, it is inferred from protein.stable.ID, but it may not be protein actually

used.genes.human <- unique(ensembl.9606.db[which(ensembl.9606.db[,2] %in% to.use.human.orig), 'Gene.name'])
res.genes.human <- unique(intersect(used.genes.human, 
	entrez.9606.db[which(entrez.9606.db$type_of_gene == "protein-coding"), "Symbol_from_nomenclature_authority"]))
# This get to have 18910 genes, aligned with previous result

# CCDC28A-AS1 is an confusing one, when it exists in *.dummy, but not one protein-coding gene
res.genes.human.dummy <- unique(intersect(used.genes.human, entrez.9606.db$Symbol_from_nomenclature_authority))


# conclusion
# Treat protein.stable.ID as protein, not regarding NCBI protein-coding tags
# Finally, get 18990 genes (inferred from protein)
#    | collected by pairs.human.db




#### Mouse part
# ensembl.10090.db get 21291 proteins (inferred from protein.stable.ID)

mouse.orig.protein <- unique(c(string.10090.readin.db[,1], string.10090.readin.db[,2]))

mouse.split.protein <- lapply(mouse.orig.protein, function(x) {strsplit(x, split=".", fixed=TRUE) })
mouse.split.protein <- unlist(mouse.split.protein)  
mouse.split.protein <- mouse.split.protein[2* (1:(length(mouse.split.protein)/2))]
mouse.split.protein <- unique(mouse.split.protein)


to.use.mouse.orig <- unique(intersect(mouse.split.protein, ensembl.10090.db[,2]))
length(to.use.mouse.orig) # 21291, it is inferred from protein.stable.ID, but it may not be protein actually

used.genes.mouse <- unique(ensembl.10090.db[which(ensembl.10090.db[,2] %in% to.use.mouse.orig), 'Gene.name'])
res.genes.mouse <- unique(intersect(used.genes.mouse, 
	entrez.10090.db[which(entrez.10090.db$type_of_gene == "protein-coding"), "Symbol_from_nomenclature_authority"]))
# This get to have 20487 genes, not aligned with previous result

res.genes.mouse.dummy <- unique(intersect(used.genes.mouse, entrez.10090.db$Symbol_from_nomenclature_authority))
# 20941, a little more than previous result

# It is different from previous result, as 
# NAs are removed in previous analysis for one gene partner NA to remove the involved gene pairs


# conclusion
# Treat protein.stable.ID as protein, not regarding NCBI protein-coding tags
# Finally, get 18990 genes (inferred from protein)
#    | collected by pairs.mouse.db


