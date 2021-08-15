# STRING database explore
# STRING use Ensembl IDs


# entrez db need to be specified in the beginning
entrez.10090.db

# ref db: Ensembl database
#
ensembl.10090.db <- read.csv("../STRING/10090-csv-mgi-mart_export.txt", stringsAsFactors = FALSE)
ensembl.10090.db <- ensembl.10090.db[which(ensembl.10090.db$Protein.stable.ID != ""), ]
# not use this line anymore at 2020.12.20
# ensembl.10090.db <- ensembl.10090.db[which(ensembl.10090.db$MGI.ID != ""), ]



### ensembl db repack-up with GchR37 & GchR38 and manual check all unmatched rest


# string PPI database
string.10090.readin.db <- read.table("../STRING/10090.protein.links.full.v11.0-20191023.txt",
					header = TRUE,
					sep = " ",
					quote = "",
					comment.char = "",
					stringsAsFactors = FALSE)

inter.CUT.SCORE <- 399
string.10090.raw.db <- string.10090.readin.db[which(string.10090.readin.db$combined_score > inter.CUT.SCORE), ]
# 773760 rows, at 2019.11.09

# to replace OLD-method
# step 1: match to ensembl
# use gene.name to directly match
# string proteins factor
pro1.factor <- levels(factor(string.10090.raw.db$protein1))
pro2.factor <- levels(factor(string.10090.raw.db$protein2))
# function
# this function is to shorten the code. make sure var: prog.bar.j exists
pro.ensembl.ID.len <- 24  # nchar(), here to show but of no direct usage
ProFactorMapGeneName <- function(x) {
	pro.ID <- substr(x, 7, 24)
	ind.match <- which(pro.ID == ensembl.10090.db$Protein.stable.ID)
	gene.name <- NA
	if (length(ind.match) > 1) {  # multi-matches, check if gene.names are the same
		equal.match <- TRUE
		for (i in 2:length(ind.match)) {
			equal.match <- equal.match && (ensembl.10090.db$Gene.name[ind.match[i-1]] == ensembl.10090.db$Gene.name[ind.match[i]])
		}
		if (!equal.match)
			warning(paste0("Multi-matches with different Gene.name in ensembl.10090.db!", ensembl.10090.db$Gene.name[ind.match[1]]))
	}
	if (length(ind.match) > 0) {  # find
		gene.name <- ensembl.10090.db$Gene.name[ind.match[1]]  # use only match[1]
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
	ind.match <- which(x == entrez.10090.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	gene.fullname <- NA
	gene.syno <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches!")
	if (length(ind.match) > 0) {
		gene.ID <- entrez.10090.db$GeneID[ind.match[1]]  # use only match[1]
		gene.fullname <- entrez.10090.db$Full_name_from_nomenclature_authority[ind.match[1]]
		gene.syno <- entrez.10090.db$Synonyms[ind.match[1]]
	}
	prog.bar.j$tick()
	c(gene.ID, gene.fullname, gene.syno)
})
# protein 2
prog.bar.j <- progress::progress_bar$new(total = length(pro2.factor))
prog.bar.j$tick(0)
pro2.factor.mapped.res <- sapply(pro2.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.10090.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	gene.fullname <- NA
	gene.syno <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches!")
	if (length(ind.match) > 0) {
		gene.ID <- entrez.10090.db$GeneID[ind.match[1]]  # use only match[1]
		gene.fullname <- entrez.10090.db$Full_name_from_nomenclature_authority[ind.match[1]]
		gene.syno <- entrez.10090.db$Synonyms[ind.match[1]]
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
string.pro1.ID <- pro1.factor.mapped.ID[match(string.10090.raw.db$protein1, pro1.factor)]
string.pro2.ID <- pro2.factor.mapped.ID[match(string.10090.raw.db$protein2, pro2.factor)]
# string db matching to Symbols
string.pro1.genename <- pro1.factor.mapped.Gene.name[match(string.10090.raw.db$protein1, pro1.factor)]
string.pro2.genename <- pro2.factor.mapped.Gene.name[match(string.10090.raw.db$protein2, pro2.factor)]
# string db matching to FullName
string.pro1.fullname <- pro1.factor.mapped.fullname[match(string.10090.raw.db$protein1, pro1.factor)]
string.pro2.fullname <- pro2.factor.mapped.fullname[match(string.10090.raw.db$protein2, pro2.factor)]
# string db matching to synonyms
string.pro1.syno <- pro1.factor.mapped.syno[match(string.10090.raw.db$protein1, pro1.factor)]
string.pro2.syno <- pro2.factor.mapped.syno[match(string.10090.raw.db$protein2, pro2.factor)]


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
interaction.list.string[, colnames.select.list] <- string.10090.raw.db[, colnames.select.list]
interaction.list.string <- interaction.list.string[which(!is.na(interaction.list.string$inter.GeneID.A)), ]
interaction.list.string <- interaction.list.string[which(!is.na(interaction.list.string$inter.GeneID.B)), ]
# whether doing unique(), depending on how many lines are the result, if it got too many, it will take a long long time
# interaction.list.string <- unique(interaction.list.string)
# rows 722440, at 2019.11.15, as mt-DNAs included
# 
# rows 720723, at 2019.11.09, after unique()
# 
# rows 11195499, at 2019.10.21
# select score over ???
# select 900, at 2019.10.21
#inter.CUT.SCORE <- 700  # defined in STRING article, as > 0.7 is of high confidence
#interaction.list.string <- interaction.list.string[which(interaction.list.string[,"combined_score"] > inter.CUT.SCORE), ]
#
# to speed up, and lower the algorithm cost, use inter.CUT.SCORE in the readin stage. See in the rather top lines.
#
# rows 11944806, at 2020.12.20
# after remove NAs, 11790260 rows remained!




# [added at 2020.12.20] from human scripts
## FOR generating actions database
#
# actions.txt is added into the database
actions.10090.readin.db <- read.table("../STRING/10090.protein.actions.v11.0.txt",
					header = TRUE,
					sep = "\t",
					quote = "",
					comment.char = "",
					stringsAsFactors = FALSE)
# 4850272 rows, at 2020.12.20
#
# ------ use the same steps like the string.raw.db ------
# step 1: match to ensembl
# use gene.name to directly match
# string proteins factor
it.a.factor <- levels(factor(actions.10090.readin.db$item_id_a))
it.b.factor <- levels(factor(actions.10090.readin.db$item_id_b))
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
	ind.match <- which(x == entrez.10090.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches! with GeneName: ", x)
	if (length(ind.match) > 0) {
		gene.ID <- entrez.10090.db$GeneID[ind.match[1]]  # use only match[1]
	}
	prog.bar.ja$tick()
	gene.ID
})
# protein 2
prog.bar.jb <- progress::progress_bar$new(total = length(it.b.factor))
prog.bar.jb$tick(0)
it.b.factor.mapped.res <- sapply(it.b.factor.mapped.Gene.name, USE.NAMES = FALSE, function(x) {
	ind.match <- which(x == entrez.10090.db$Symbol_from_nomenclature_authority)
	gene.ID <- NA
	if (length(ind.match) > 1)
		warning("unexpected multi-matches!")
	if (length(ind.match) > 0) {
		gene.ID <- entrez.10090.db$GeneID[ind.match[1]]  # use only match[1]
	}
	prog.bar.jb$tick()
	gene.ID
})
#
# actions.db matching to IDs
match.it.a.ID <- it.a.factor.mapped.res[match(actions.10090.readin.db$item_id_a, it.a.factor)]
match.it.b.ID <- it.b.factor.mapped.res[match(actions.10090.readin.db$item_id_b, it.b.factor)]
# actions.db matching to Symbols
match.it.a.genename <- it.a.factor.mapped.Gene.name[match(actions.10090.readin.db$item_id_a, it.a.factor)]
match.it.b.genename <- it.b.factor.mapped.Gene.name[match(actions.10090.readin.db$item_id_b, it.b.factor)]
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
interaction.attached.actions[, colnames.select.remained] <- actions.10090.readin.db[, colnames.select.remained]
interaction.attached.actions <- interaction.attached.actions[which(!is.na(interaction.attached.actions$inter.GeneID.A)), ]
interaction.attached.actions <- interaction.attached.actions[which(!is.na(interaction.attached.actions$inter.GeneID.B)), ]
# 4850272 rows -> 4763978 rows, at 2020.12.20
#



#
# below
#
# NOT RUNNING, just as test
#
#
# using following to slim string.db
# GO-list
# Ligand
go.10090.ligand.db
# Receptor
go.10090.receptor.db
# ECM
go.10090.ecm.db

# Mesh-HGNC-list
# Ligand
hgnc.10090.ligand.cut.db
# Receptor
hgnc.10090.receptor.cut.db
# ECM
hgnc.10090.ecm.cut.db

# GeneID list (LRE)
go.GeneID.list <- c(go.10090.ligand.db$GeneID, go.10090.receptor.db$GeneID, go.10090.ecm.db$GeneID)
hgnc.GeneID.list <- c(hgnc.10090.ligand.cut.db$GeneID, hgnc.10090.receptor.cut.db$GeneID, hgnc.10090.ecm.cut.db$GeneID)
#
LRE.GeneID.list <- c(go.GeneID.list, hgnc.GeneID.list)
LRE.GeneID.list <- unique(LRE.GeneID.list)
#
# slimming string.db
interaction.list.string <- interaction.list.string[which(!is.na(match(interaction.list.string$inter.GeneID.A, LRE.GeneID.list))),]
interaction.list.string <- interaction.list.string[which(!is.na(match(interaction.list.string$inter.GeneID.B, LRE.GeneID.list))),]
interaction.list.string <- unique(interaction.list.string)
# rows 61622, at 2019.10.21
# use > 700, get rows 101739, at 2019.10.22
# use > 700, get rows 101743, at 2019.11.09








