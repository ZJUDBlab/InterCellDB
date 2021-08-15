
# entrez db need to be specified in the beginning
entrez.9606.db

# ref db: Ensembl database
#
ensembl.9606.db <- read.csv("../STRING/ensembl-9606-20191115.csv", stringsAsFactors = FALSE)
ensembl.9606.db <- ensembl.9606.db[which(ensembl.9606.db$Protein.stable.ID != ""), ]
#### consider not used this line!!!
#ensembl.9606.db <- ensembl.9606.db[which(!is.na(ensembl.9606.db$NCBI.gene.ID)), ]
# at 2019.11.15, 109689 rows

### ensembl db repack-up with GchR37 & GchR38 and manual check all unmatched rest
# 131631 rows, at 2020.12.20

# string PPI database
string.9606.readin.db <- read.table("../STRING/9606.protein.links.full.v11.0-20191115.txt",
					header = TRUE,
					sep = " ",
					quote = "",
					comment.char = "",
					stringsAsFactors = FALSE)
inter.CUT.SCORE <- 700
string.9606.raw.db <- string.9606.readin.db[which(string.9606.readin.db$combined_score > inter.CUT.SCORE), ]
# 839522 rows, at 2019.11.15
# 
# select which(experiments.score != 0)
string.9606.exp.only.db <- string.9606.readin.db[which(string.9606.readin.db$experiments != 0), ]

# step 1: match to ensembl
# use gene.name to directly match
# string proteins factor
pro1.factor <- levels(factor(string.9606.raw.db$protein1))
pro2.factor <- levels(factor(string.9606.raw.db$protein2))
# function
# this function is to shorten the code. make sure var: prog.bar.j exists
pro.ensembl.ID.len <- 20  # nchar(), here to show but of no direct usage
ProFactorMapGeneName <- function(x) {
	pro.ID <- substr(x, 6, 20)
	ind.match <- which(pro.ID == ensembl.9606.db$Protein.stable.ID)
	gene.name <- NA
	if (length(ind.match) > 1) {  # multi-matches, check if gene.names are the same
		equal.match <- TRUE
		for (i in 2:length(ind.match)) {
			equal.match <- equal.match && (ensembl.9606.db$Gene.name[ind.match[i-1]] == ensembl.9606.db$Gene.name[ind.match[i]])
		}
		if (!equal.match) {
			warning(paste0("Multi-matches with different Gene.name in ensembl.9606.db!", ensembl.9606.db$Gene.name[ind.match[1]]))
		}
	}
	if (length(ind.match) > 0) {  # find
		gene.name <- ensembl.9606.db$Gene.name[ind.match[1]]  # use only match[1]
	}
	prog.bar.j$tick()
	gene.name
}
# protein 1
prog.bar.j <- progress::progress_bar$new(total = length(pro1.factor))
prog.bar.j$tick(0)
pro1.factor.mapped.Gene.name <- sapply(pro1.factor, USE.NAMES = FALSE, ProFactorMapGeneName)
# protein 2
prog.bar.j <- progress::progress_bar$new(total = length(pro2.factor))
prog.bar.j$tick(0)
pro2.factor.mapped.Gene.name <- sapply(pro2.factor, USE.NAMES = FALSE, ProFactorMapGeneName)
#
#
# step 2: match to entrez IDs & synonyms
# protein 1
prog.bar.j <- progress::progress_bar$new(total = length(pro1.factor))
prog.bar.j$tick(0)
pro1.factor.mapped.res <- sapply(pro1.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.9606.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	gene.fullname <- NA
	gene.syno <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches! with GeneName: ", x)
	if (length(ind.match) > 0) {
		gene.ID <- entrez.9606.db$GeneID[ind.match[1]]  # use only match[1]
		gene.fullname <- entrez.9606.db$Full_name_from_nomenclature_authority[ind.match[1]]
		gene.syno <- entrez.9606.db$Synonyms[ind.match[1]]
	}
	prog.bar.j$tick()
	c(gene.ID, gene.fullname, gene.syno)
})
# protein 2
prog.bar.j <- progress::progress_bar$new(total = length(pro2.factor))
prog.bar.j$tick(0)
pro2.factor.mapped.res <- sapply(pro2.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.9606.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	gene.fullname <- NA
	gene.syno <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches!")
	if (length(ind.match) > 0) {
		gene.ID <- entrez.9606.db$GeneID[ind.match[1]]  # use only match[1]
		gene.fullname <- entrez.9606.db$Full_name_from_nomenclature_authority[ind.match[1]]
		gene.syno <- entrez.9606.db$Synonyms[ind.match[1]]
	}
	prog.bar.j$tick()
	c(gene.ID, gene.fullname, gene.syno)
})
#
# ID
pro1.factor.mapped.ID <- t(pro1.factor.mapped.res)[,1]
pro2.factor.mapped.ID <- t(pro2.factor.mapped.res)[,1]
# fullname
pro1.factor.mapped.fullname <- t(pro1.factor.mapped.res)[,2]
pro2.factor.mapped.fullname <- t(pro2.factor.mapped.res)[,2]
# synonyms
pro1.factor.mapped.syno <- t(pro1.factor.mapped.res)[,3]
pro2.factor.mapped.syno <- t(pro2.factor.mapped.res)[,3]
#
# string db matching to IDs
# as string db contains too many rows, use grep() match all approved IDs
string.pro1.ID <- pro1.factor.mapped.ID[match(string.9606.raw.db$protein1, pro1.factor)]
string.pro2.ID <- pro2.factor.mapped.ID[match(string.9606.raw.db$protein2, pro2.factor)]
# string db matching to Symbols
string.pro1.genename <- pro1.factor.mapped.Gene.name[match(string.9606.raw.db$protein1, pro1.factor)]
string.pro2.genename <- pro2.factor.mapped.Gene.name[match(string.9606.raw.db$protein2, pro2.factor)]
# string db matching to FullName
string.pro1.fullname <- pro1.factor.mapped.fullname[match(string.9606.raw.db$protein1, pro1.factor)]
string.pro2.fullname <- pro2.factor.mapped.fullname[match(string.9606.raw.db$protein2, pro2.factor)]
# string db matching to synonyms
string.pro1.syno <- pro1.factor.mapped.syno[match(string.9606.raw.db$protein1, pro1.factor)]
string.pro2.syno <- pro2.factor.mapped.syno[match(string.9606.raw.db$protein2, pro2.factor)]


#
# construct interaction pairs
interaction.list.string <- data.frame(
	inter.GeneID.A = string.pro1.ID,
	inter.GeneID.B = string.pro2.ID,
	inter.GeneName.A = string.pro1.genename,
	inter.GeneName.B = string.pro2.genename,
	inter.FullName.A = string.pro1.fullname,
	inter.FullName.B = string.pro2.fullname,
	inter.Synonyms.A = string.pro1.syno,
	inter.Synonyms.B = string.pro2.syno,
	stringsAsFactors = FALSE)
colnames.select.list <- c("neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
						  "homology", "coexpression", "coexpression_transferred", "experiments",
						  "experiments_transferred", "database", "database_transferred", 
						  "textmining", "textmining_transferred", "combined_score")
interaction.list.string[, colnames.select.list] <- string.9606.raw.db[, colnames.select.list]
interaction.list.string <- interaction.list.string[which(!is.na(interaction.list.string$inter.GeneID.A)), ]
interaction.list.string <- interaction.list.string[which(!is.na(interaction.list.string$inter.GeneID.B)), ]
# 779596 rows, at 2019.11.15
# generate experiments-only (no score limit)
# 417192 rows, at 2019.11.30
#
# 11759454 rows, with NA removed, there remains 11096258 rows at 2020.11.28
# 11759454 rows -> 11560956 rows at 2020.12.20, using new ensembl databse



# [new addition at 2019.11.29]
# actions.txt is added into the database
actions.9606.readin.db <- read.table("../STRING/9606.protein.actions.v11.0-20191122.txt",
					header = TRUE,
					sep = "\t",
					quote = "",
					comment.char = "",
					stringsAsFactors = FALSE)
# 3470906 rows, at 2019.11.29
#
# ------ use the same steps like the string.raw.db ------
#
# step 1: match to ensembl
# use gene.name to directly match
# string proteins factor
it.a.factor <- levels(factor(actions.9606.readin.db$item_id_a))
it.b.factor <- levels(factor(actions.9606.readin.db$item_id_b))
#
# [function imported] 
# importFrom the above ProFactorMapGeneName <- function(x) {}
#
# protein 1
prog.bar.j <- progress::progress_bar$new(total = length(it.a.factor))
prog.bar.j$tick(0)
it.a.factor.mapped.Gene.name <- sapply(it.a.factor, USE.NAMES = FALSE, ProFactorMapGeneName)
# protein 2
prog.bar.j <- progress::progress_bar$new(total = length(it.b.factor))
prog.bar.j$tick(0)
it.b.factor.mapped.Gene.name <- sapply(it.b.factor, USE.NAMES = FALSE, ProFactorMapGeneName)
#
#
# step 2: match to entrez IDs & synonyms
# protein 1
prog.bar.ja <- progress::progress_bar$new(total = length(it.a.factor))
prog.bar.ja$tick(0)
it.a.factor.mapped.res <- sapply(it.a.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.9606.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches! with GeneName: ", x)
	if (length(ind.match) > 0) {
		gene.ID <- entrez.9606.db$GeneID[ind.match[1]]  # use only match[1]
	}
	prog.bar.ja$tick()
	gene.ID
})
# protein 2
prog.bar.jb <- progress::progress_bar$new(total = length(it.b.factor))
prog.bar.jb$tick(0)
it.b.factor.mapped.res <- sapply(it.b.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.9606.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches!")
	if (length(ind.match) > 0) {
		gene.ID <- entrez.9606.db$GeneID[ind.match[1]]  # use only match[1]
	}
	prog.bar.jb$tick()
	gene.ID
})
#
# actions.db matching to IDs
match.it.a.ID <- it.a.factor.mapped.res[match(actions.9606.readin.db$item_id_a, it.a.factor)]
match.it.b.ID <- it.b.factor.mapped.res[match(actions.9606.readin.db$item_id_b, it.b.factor)]
# actions.db matching to Symbols
match.it.a.genename <- it.a.factor.mapped.Gene.name[match(actions.9606.readin.db$item_id_a, it.a.factor)]
match.it.b.genename <- it.b.factor.mapped.Gene.name[match(actions.9606.readin.db$item_id_b, it.b.factor)]
#
# construct actions-annotated interaction pairs
interaction.attached.actions <- data.frame(
	inter.GeneID.A = match.it.a.ID,
	inter.GeneID.B = match.it.b.ID,
	inter.GeneName.A = match.it.a.genename,
	inter.GeneName.B = match.it.b.genename,
	stringsAsFactors = FALSE)
colnames.select.remained <- c("mode", "action", "is_directional", "a_is_acting", "score")
# add rest cols
interaction.attached.actions[, colnames.select.remained] <- actions.9606.readin.db[, colnames.select.remained]
interaction.attached.actions <- interaction.attached.actions[which(!is.na(interaction.attached.actions$inter.GeneID.A)), ]
interaction.attached.actions <- interaction.attached.actions[which(!is.na(interaction.attached.actions$inter.GeneID.B)), ]
# 3139858 rows, at 2019.11.29
# 3470906 rows -> 3409770 rows, at 2020.12.20


#
# as far as all we need is the LRE-related database
# Do data-slim based on LRE-annotation
#
# use ligand.anno.m.db, recepetor.anno.m.db, ecm.anno.m.db to construct LRE.anno.m.db
LRE.anno.m.db  # 4419 rows, at 2019.11.30
#
actions.9606.db # <- interaction.attached.actions
inds.A <- which(actions.9606.db$inter.GeneName.A %in% LRE.anno.m.db$Symbol_from_nomenclature_authority)
inds.B <- which(actions.9606.db$inter.GeneName.B %in% LRE.anno.m.db$Symbol_from_nomenclature_authority)
#
actions.9606.LRE.db <- actions.9606.db[intersect(inds.A, inds.B), ]
# 481514 rows, at 2019.11.30

# For Mouse, run at 2020.01.20
# 496206 rows



