library(dplyr)



# readin GO db (all taxonomies are included)
GO.all.db <- read.table("../For-GO/gene2go.txt",
	header = TRUE,
	sep = "\t",
	quote = "",
	comment.char = "",
	stringsAsFactors = FALSE)
colnames(GO.all.db)[1] <- "tax_id"

# ref entrez db
entrez.9606.db   # human
entrez.10090.db  # mouse


# running human GO generation process
GO.9606.raw.db <- GO.all.db[which(GO.all.db$tax_id == 9606), ]  # 323030 rows
GO.9606.pending.db <- left_join(GO.9606.raw.db, entrez.9606.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
	by = c("GeneID" = "GeneID"))
colnames(GO.9606.pending.db)[ncol(GO.9606.pending.db)] <- "Gene.name"
GO.9606.pending.db <- GO.9606.pending.db[, c(1:2, ncol(GO.9606.pending.db), 3:(ncol(GO.9606.pending.db) - 1))]
GO.9606.pending.db <- GO.9606.pending.db[which(!is.na(GO.9606.pending.db$Gene.name)), ]
# 322873 rows



# running mouse GO generation process
GO.10090.raw.db <- GO.all.db[which(GO.all.db$tax_id == 10090), ]  # 367376 rows
GO.10090.pending.db <- left_join(GO.10090.raw.db, entrez.10090.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
	by = c("GeneID" = "GeneID"))
colnames(GO.10090.pending.db)[ncol(GO.10090.pending.db)] <- "Gene.name"
GO.10090.pending.db <- GO.10090.pending.db[, c(1:2, ncol(GO.10090.pending.db), 3:(ncol(GO.10090.pending.db) - 1))]
GO.10090.pending.db <- GO.10090.pending.db[which(!is.na(GO.10090.pending.db$Gene.name)), ]
# 367110 rows

