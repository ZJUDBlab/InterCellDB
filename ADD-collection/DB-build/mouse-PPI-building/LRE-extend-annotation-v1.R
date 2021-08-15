# LRE-extend-annotation-v1.R

# step 0: collect all now generated
# PPI pairs
# known (database or experiments-based)
# ggpairs -> bbpairs when STRING added, deprecated at 2019.11.09
bbpairs


# known (database or experiments-based)
# use only db from STRING at 2019.10.29 tmp
sspairs # rows 546117

# predicted and other (from STRING db)
aapairs # rows 73198


### LRE
## GO-list
# Ligand
go.10090.ligand.db  # 3042 rows, non-unique in only cols in select.col.list
# Receptor
go.10090.receptor.db  # 505 rows, non-unique in only cols in select.col.list
# ECM
go.10090.ecm.db  # 964 rows, non-unique in only cols in select.col.list

## Mesh-HGNC-list
# Ligand
hgnc.10090.ligand.cut.db  # 328 rows
# Receptor
hgnc.10090.receptor.cut.db  # 959 rows
# ECM
hgnc.10090.ecm.cut.db  # 117 rows

# [if it is the first time]
# step 1: merge LRE list
# To achieve usage best in future, use minimally conserved data
select.col.list <- c("GeneID", "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority")
ligand.m.db <- rbind(go.10090.ligand.db[, select.col.list], hgnc.10090.ligand.cut.db[, select.col.list])
rownames(ligand.m.db) <- NULL
ligand.m.db <- unique(ligand.m.db)
# done, ligand, 1820 rows at 2019.10.26

receptor.m.db <- rbind(go.10090.receptor.db[, select.col.list], hgnc.10090.receptor.cut.db[, select.col.list])
rownames(receptor.m.db) <- NULL
receptor.m.db <- unique(receptor.m.db)
# done, receptor, 1114 rows at 2019.10.26

ecm.m.db <- rbind(go.10090.ecm.db[, select.col.list], hgnc.10090.ecm.cut.db[, select.col.list])
rownames(ecm.m.db) <- NULL
ecm.m.db <- unique(ecm.m.db)
# done, ecm, 519 rows at 2019.10.26

# step 1.adv: LRE-annotated
ligand.anno.m.db	<- readRDS("../RUNOUT/OUTs-20200130/Ligand-anno-merged-upd-20200131.rds")
receptor.anno.m.db 	<- readRDS("../RUNOUT/OUTs-20200130/Receptor-anno-merged-upd-20200131.rds")
ecm.anno.m.db 		<- readRDS("../RUNOUT/OUTs-20200130/ECM-anno-merged-upd-20200131.rds")

# check
{
	# LRE database cross check
	LRE.db <- rbind(ligand.anno.m.db, receptor.anno.m.db, ecm.anno.m.db)
	LRE.db$UID <- 1:nrow(LRE.db)
	#
	LRE.total.u.db <- unique(LRE.db)  # check if nrow() is not the same as LRE.db
	LRE.part.u.db <- DoPartUnique(LRE.db, c(1,2))
	# find those genes are included in multiple databases
	inds.dup <- setdiff(LRE.total.u.db[, "UID"], LRE.part.u.db[, "UID"])
	#
	geneID.dup <- integer()
	for (i in inds.dup) {
		this.geneid <- LRE.total.u.db[i, "GeneID"]
		inds.m <- which(this.geneid == LRE.total.u.db[, "GeneID"])
		if (length(inds.m) < 2) {
			warning(paste0("Unexpected error on GeneID: ", this.geneid))
		}
		# do unique on this
		tmp.mdb <- LRE.total.u.db[inds.m, ]
		tmp.mdb$UID <- NULL  # remove UID
		tmp.mdb <- unique(tmp.mdb)
		if (nrow(tmp.mdb) > 1) {
			geneID.dup <- c(geneID.dup, as.integer(this.geneid))
		}
	}

	check.geneid <- 114249
	ligand.anno.m.db[which(ligand.anno.m.db$GeneID == check.geneid), ]
	receptor.anno.m.db[which(receptor.anno.m.db$GeneID == check.geneid), ]
	ecm.anno.m.db[which(ecm.anno.m.db$GeneID == check.geneid), ]
	# 114249: Npnt, <Location>, <Species>
}


# step 2: recalculate the presumed score and database-exp-based score
# calculation fomula see below site: 
# https://stringdb-static.org/download/combine_subscores.py
#
RecalcScores <- function(xxpairs.score.part) 
{
	known.score.list 	<- c("database", "experiments")
	presumed.score.list <- c("neighborhood", "neighborhood_transferred", "fusion", "coexpression",
							 "coexpression_transferred", "experiments_transferred", "database_transferred")
	txtm.score.list <- c("textmining", "textmining_transferred")
	needhomo.cooccur <- c("cooccurence")
	homology.colname <- c("homology")

	kp.prior <- 0.041  # specified by STRING, [TODO] further recheck where it comes from
	ComputePriorAway <- function(score, prior) {
		if (score < prior) { score <- prior }
		score.no.prior <- (score - prior) / (1.0 - prior)
		# return
		score.no.prior
	}

	prog.bar.recalc <- progress::progress_bar$new(total = nrow(xxpairs.score.part))
	prog.bar.recalc$tick(0)
	comb.parts.scores <- apply(xxpairs.score.part, MARGIN = 1, function(x) {
		# part 1: experiments score
		score.exp.direct <- as.integer(x["experiments"])
		# part 2: database score
		score.db.direct <- as.integer(x["database"])
		# part 3: presumed score
		score.presumed.nx <- as.single(x[presumed.score.list]) / 1000
		score.presumed.nx.no.prior <- sapply(score.presumed.nx, prior = kp.prior, ComputePriorAway)
		score.presumed.res <- 1.0
		for (i in 1:length(score.presumed.nx.no.prior)) {
			score.presumed.res <- score.presumed.res * (1.0 - score.presumed.nx.no.prior[i])
		}
		# homo correction part
		score.homo <- as.single(x["homology"]) / 1000
		# cooccurence
		score.cooccur <- as.single(x["cooccurence"]) / 1000
		score.cooccur.no.prior <- ComputePriorAway(score.cooccur, kp.prior)
		score.cooccur.np.homo <- score.cooccur.no.prior * (1.0 - score.homo)
		# text mining & its transferred one
		score.txtm <- as.single(x[txtm.score.list]) / 1000
		score.txtm.no.prior <- sapply(score.txtm, prior = kp.prior, ComputePriorAway)
		score.txtm.res <- 1.0
		for (i in 1:length(score.txtm.no.prior)) {
			score.txtm.res <- score.txtm.res * (1.0 - score.txtm.no.prior[i])
		}
		score.txtm.res.homo <- (1.0 - score.txtm.res) * (1.0 - score.homo)
		# calc presumed final
		score.presumed.res.onem <- score.presumed.res * (1.0 - score.cooccur.np.homo) * (1.0 - score.txtm.res.homo)
		score.presumed.res.final <- (1.0 - score.presumed.res.onem) * (1.0 - kp.prior) + kp.prior
		score.presumed.res.final <- as.integer(1000 * score.presumed.res.final)
		# return of apply
		prog.bar.recalc$tick()
		score.res3 <- if (score.presumed.res.final > 1000 * kp.prior) score.presumed.res.final else 0
		score.res3 <- as.integer(score.res3)
		c(score.exp.direct, score.db.direct, score.res3)
	})
	# return
	comb.parts.scores <- t(comb.parts.scores)
}

recalc.sel.cols <- c("neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
					 "homology", "coexpression", "coexpression_transferred", "experiments",
					 "experiments_transferred", "database", "database_transferred", 
					 "textmining", "textmining_transferred", "combined_score")

aapairs.score.part <- as.matrix(aapairs[, recalc.sel.cols])
aapairs.recalc.score.part <- RecalcScores(aapairs.score.part)
aapairs$inter.Experiments.Score <- aapairs.recalc.score.part[, 1]
aapairs$inter.Database.Score <- aapairs.recalc.score.part[, 2]
aapairs$inter.Predicted.Score <- aapairs.recalc.score.part[, 3]
aapairs$inter.Combined.Score <- aapairs[, "combined_score"]
aapairs.recalc <- aapairs[, c(1:8, 23:26)]

sspairs.score.part <- as.matrix(sspairs[, recalc.sel.cols])
sspairs.recalc.score.part <- RecalcScores(sspairs.score.part)
sspairs$inter.Experiments.Score <- sspairs.recalc.score.part[, 1]
sspairs$inter.Database.Score <- sspairs.recalc.score.part[, 2]
sspairs$inter.Predicted.Score <- sspairs.recalc.score.part[, 3]
sspairs$inter.Combined.Score <- sspairs[, "combined_score"]
sspairs.recalc <- sspairs[, c(1:8, 23:26)]

# alert at 2020.11.28, STRING get gene pairs has different scores
# e.g. AK2~ADSL, 2 records are presented, 
# 	1. AK2~ADSL Combined_score = 289, database = 0
#	2. AK2~ADSL Combined_score = 967, database = 900


# step 3: construct LRE database
#
# @param xxpairs. data.frame. format must be like bbpairs or aapairs
# @param ligand.m.db		format like the above
# @param receptor.m.db		format like the above
# @param ecm.m.db			format like the above
# @param ifanno		logic. if TRUE, add annotation infos in the return value
LREMapPairs <- function(xxpairs, ligand.m.db, receptor.m.db, ecm.m.db, ifanno = FALSE) 
{
 # use 3 binary code
 # [binary]	2		1			0
 #			ligand	receptor	ecm
 #
 #
 LRE.bin.list <- c(
	NA, # 0
	"ECM", 	# 1
	"Receptor", # 2
	"Receptor, ECM", # 3
	"Ligand", # 4
	"Ligand, ECM", # 5
	"Ligand, Receptor", # 6
	"Ligand, Receptor, ECM" # 7
 )
 # ifanno = TRUE, use informations got from thes columns
 annotation.infos <- c("Location", "Type", "Species")
 	#
	### xxpairs LRE apart
	bb.inter.ID.A <- levels(factor(xxpairs$inter.GeneID.A))
	bb.inter.ID.B <- levels(factor(xxpairs$inter.GeneID.B))
	# inter A LRE mapping
	prog.bar.LRE.A.packp <- progress::progress_bar$new(total = length(bb.inter.ID.A))
	prog.bar.LRE.A.packp$tick(0)
	bb.inter.A.LRE.pack <- sapply(bb.inter.ID.A, USE.NAMES = FALSE, function(x) {
		ind.m.ligand	<- match(x, ligand.m.db$GeneID)
		ind.m.receptor 	<- match(x, receptor.m.db$GeneID)
		ind.m.ecm		<- match(x, ecm.m.db$GeneID)
		ligand.ifnot 	<- !is.na(ind.m.ligand)
		receptor.ifnot 	<- !is.na(ind.m.receptor)
		ecm.ifnot 		<- !is.na(ind.m.ecm)
		calc.lre.val <- (if (ligand.ifnot) 4 else 0) + (if (receptor.ifnot) 2 else 0) + (if (ecm.ifnot) 1 else 0)
		prog.bar.LRE.A.packp$tick()
		c(calc.lre.val, ind.m.ligand, ind.m.receptor, ind.m.ecm)
	})
	bb.inter.A.LRE.pack <- t(bb.inter.A.LRE.pack)  # why t(), as it is add-up by cols, so the same meaning ones is in the same row(if not tranpose)
	# inter B LRE mapping
	prog.bar.LRE.B.packp <- progress::progress_bar$new(total = length(bb.inter.ID.B))
	prog.bar.LRE.B.packp$tick(0)
	bb.inter.B.LRE.pack <- sapply(bb.inter.ID.B, USE.NAMES = FALSE, function(x) {
		ind.m.ligand	<- match(x, ligand.m.db$GeneID)
		ind.m.receptor 	<- match(x, receptor.m.db$GeneID)
		ind.m.ecm		<- match(x, ecm.m.db$GeneID)
		ligand.ifnot 	<- !is.na(ind.m.ligand)
		receptor.ifnot 	<- !is.na(ind.m.receptor)
		ecm.ifnot 		<- !is.na(ind.m.ecm)
		calc.lre.val <- (if (ligand.ifnot) 4 else 0) + (if (receptor.ifnot) 2 else 0) + (if (ecm.ifnot) 1 else 0)
		prog.bar.LRE.B.packp$tick()
		c(calc.lre.val, ind.m.ligand, ind.m.receptor, ind.m.ecm)
	})
	bb.inter.B.LRE.pack <- t(bb.inter.B.LRE.pack)
	# mapped to xxpairs
	xxpairs.mapped.to.bb.inter.A <- match(xxpairs$inter.GeneID.A, bb.inter.ID.A)
	xxpairs.mapped.to.bb.inter.B <- match(xxpairs$inter.GeneID.B, bb.inter.ID.B)
	# part1: LRE val
	bb.inter.A.LRE.val.mapping <- bb.inter.A.LRE.pack[xxpairs.mapped.to.bb.inter.A, 1]
	bb.inter.B.LRE.val.mapping <- bb.inter.B.LRE.pack[xxpairs.mapped.to.bb.inter.B, 1]
	# part2(OPTIONAL): other info mapping and check (Location, Type(s), Species)
	if (ifanno) {
	LREMappingExtractOtherInfos <- function(xxpairs.mapped.to.bb.inter.X, bb.inter.X.LRE.pack, ligand.m.db, receptor.m.db, ecm.m.db) {
		prog.bar.LRE.X.infoget <- progress::progress_bar$new(total = length(xxpairs.mapped.to.bb.inter.X))
		prog.bar.LRE.X.infoget$tick(0)
		bb.inter.X.LRE.info.mapping <- sapply(xxpairs.mapped.to.bb.inter.X, pack = bb.inter.X.LRE.pack, 
			ligand.m.db = ligand.m.db, receptor.m.db = receptor.m.db, ecm.m.db = ecm.m.db, 
		function(x, pack, ligand.m.db, receptor.m.db, ecm.m.db) {
			ind.m.ligand 	<- pack[x, 2]
			ind.m.receptor 	<- pack[x, 3]
			ind.m.ecm 		<- pack[x, 4]
			info.getfrom.ligand 	<- if (!is.na(ind.m.ligand))	ligand.m.db[ind.m.ligand, annotation.infos]
										else c("", "", "")
			info.getfrom.receptor	<- if (!is.na(ind.m.receptor)) 	receptor.m.db[ind.m.receptor, annotation.infos]
										else c("", "", "")
			info.getfrom.ecm		<- if (!is.na(ind.m.ecm))		ecm.m.db[ind.m.ecm, annotation.infos]
										else c("", "", "")
			# check if all sourced infos are the same
			check.lr <- TRUE
			check.le <- TRUE
			check.re <- TRUE
			if (!is.na(ind.m.ligand) && !is.na(ind.m.receptor)) {
				check.lr <- if (sum(info.getfrom.ligand == info.getfrom.receptor) == length(info.getfrom.ligand)) TRUE
							else FALSE
			}
			if (!is.na(ind.m.ligand) && !is.na(ind.m.ecm)) {
				check.le <- if (sum(info.getfrom.ligand == info.getfrom.ecm) == length(info.getfrom.ecm)) TRUE
							else FALSE
			}
			if (!is.na(ind.m.receptor) && !is.na(ind.m.ecm)) {
				check.re <- if (sum(info.getfrom.receptor == info.getfrom.ecm) == length(info.getfrom.receptor)) TRUE
							else FALSE
			}
			if (!(check.lr && check.le && check.re)) {  # some wrong in db generation
				warning(paste0("Cannot match infos from different sources, these are Ligand: ", ind.m.ligand, ", Receptor: ", ind.m.receptor, ", ECM: ", ind.m.ecm, "."))
			}
			info.res <- if (!is.na(ind.m.ligand)) info.getfrom.ligand
						else {
							if (!is.na(ind.m.receptor)) info.getfrom.receptor
							else {
								info.getfrom.ecm  # whether is TRUE doesn't matter, as c("","","") specified as default value
							}
						}
			prog.bar.LRE.X.infoget$tick()
			info.res <- as.character(info.res)
		})
		bb.inter.X.LRE.info.mapping <- t(bb.inter.X.LRE.info.mapping)
	}
		# with 3 cols, a serial of "Location", "Type", "Species"
		bb.inter.A.LRE.otherinfo.mapping <- LREMappingExtractOtherInfos(xxpairs.mapped.to.bb.inter.A, bb.inter.A.LRE.pack, ligand.m.db, receptor.m.db, ecm.m.db)
		bb.inter.B.LRE.otherinfo.mapping <- LREMappingExtractOtherInfos(xxpairs.mapped.to.bb.inter.B, bb.inter.B.LRE.pack, ligand.m.db, receptor.m.db, ecm.m.db)
		# append to xxpairs	
		xxpairs.ext <- cbind(xxpairs, 
						inter.LREval.A = bb.inter.A.LRE.val.mapping, 
						inter.LREval.B = bb.inter.B.LRE.val.mapping,
						inter.Location.A = bb.inter.A.LRE.otherinfo.mapping[, 1],
						inter.Location.B = bb.inter.B.LRE.otherinfo.mapping[, 1],
						inter.Type.A = bb.inter.A.LRE.otherinfo.mapping[, 2],
						inter.Type.B = bb.inter.B.LRE.otherinfo.mapping[, 2],
						inter.Species.A = bb.inter.A.LRE.otherinfo.mapping[, 3],
						inter.Species.B = bb.inter.B.LRE.otherinfo.mapping[, 3],
						stringsAsFactors = FALSE
						)
	} else {
		# append to xxpairs	
		xxpairs.ext <- cbind(xxpairs, 
						inter.LREval.A = bb.inter.A.LRE.val.mapping, 
						inter.LREval.B = bb.inter.B.LRE.val.mapping,
						stringsAsFactors = FALSE
						)
	}
	# remove 0s
	ind.ext.NAs <- union(which(xxpairs.ext$inter.LREval.A == 0), which(xxpairs.ext$inter.LREval.B == 0))
	if(length(ind.ext.NAs) > 0)
		xxpairs.ext <- xxpairs.ext[-ind.ext.NAs, ]
	# generating final db
	xxpairs.final <- xxpairs.ext
	xxpairs.final$inter.LREtype.A <- LRE.bin.list[xxpairs.final$inter.LREval.A + 1]  # step over NA
	xxpairs.final$inter.LREtype.B <- LRE.bin.list[xxpairs.final$inter.LREval.B + 1]  # step over NA
	# return
	#xxpairs.final <- unique(xxpairs.final) # don't do
	xxpairs.final
}

# this function is to do unique() only in several columes
# e.g.
#	122 134	200
#	122 134 300
# above 2 rows are the same if only unique([, c(1,2)])
DoPartUnique <- function(xxpairs, cols.select=c(1:2)) 
{
	if (sum(cols.select %in% c(1:ncol(xxpairs))) != length(cols.select)) {
		stop("Columns selected are undefined! Please check again!")
	}
	# rownames(xxpairs) <- NULL
	tmp.uni <- xxpairs[, cols.select]
	rownames(tmp.uni) <- NULL
	tmp.uni <- unique(tmp.uni)
	xxpairs <- xxpairs[as.integer(rownames(tmp.uni)),]
	xxpairs
}

# this function is to re-order the interID, and make the less one in the *.A
# e.g.
# inter.GeneID.A inter.GeneID.B
#	10				8
# after this function
# inter.GeneID.A inter.GeneID.B
#	8				10
#
# @param ind.colname.end.dual when xxpairs ends column like "xx.A xx.B"
FastAlignPairs <- function(xxpairs, ind.colname.end.dual) 
{
	inds.rev <- which(as.integer(xxpairs[, 1]) >= as.integer(xxpairs[, 2]))
	xxpairs.result <- NULL
	if (length(inds.rev) > 0) {
		xxpairs.rev <- xxpairs[inds.rev, ]
		xxpairs.conv <- xxpairs[-inds.rev, ]
		if (ncol(xxpairs.rev) <= ind.colname.end.dual) {
			if (ncol(xxpairs.rev) < ind.colname.end.dual) {
				stop("Given database has less cols than given parameter: ", ind.colname.end.dual)
			}
			xxpairs.rev.rem <- xxpairs.rev[, c(ReverseOddEvenCols(ind.colname.end.dual))]
		} else {
			xxpairs.rev.rem <- xxpairs.rev[, c(ReverseOddEvenCols(ind.colname.end.dual), (ind.colname.end.dual+1):ncol(xxpairs.rev))]
		}
		colnames(xxpairs.rev.rem) <- colnames(xxpairs.conv)
		xxpairs.result <- rbind(xxpairs.conv, xxpairs.rev.rem)
	} else {
		# do nothing and return	
		xxpairs.result <- xxpairs
	}
	# return
	xxpairs.result
}

# this function is to remove duplicate interaction pairs whether it is reversed same
# ATTENTION: this function must be used after DoPartUnique()
# e.g. 
#	122 134
#	134 122
# above 2 rows are the same, and will be removed by this function
RemoveDupPairs <- function(xxpairs) 
{
	prog.bar.rdp.1match <- progress::progress_bar$new(total = nrow(xxpairs))
	prog.bar.rdp.1match$tick(0)
	xxpairs.rev.match <- apply(xxpairs, MARGIN = 1, pairs.rev = xxpairs, function(x, pairs.rev) {
		ind.rev.m1 <- which(x[1] == pairs.rev[, 2])
		ind.rev.m2 <- which(x[2] == pairs.rev[, 1])
		ind.revm <- intersect(ind.rev.m1, ind.rev.m2)
		ifsame.log <- if (x[1] == x[2]) TRUE else FALSE
		ind.rev.res <- NA
		if (length(ind.revm) > 0 && !ifsame.log) {
			ind.rev.res <- ind.revm
		}
		prog.bar.rdp.1match$tick()
		ind.rev.res
	})
	xxpairs.rev.match <- as.integer(unlist(xxpairs.rev.match))  # use this line to check if the match part is only one-by-one, i.e. xxpairs have been unique()-ed
	xxpairs.conv.rev.raw <- data.frame(
								revID = xxpairs.rev.match, 
								convID = 1:nrow(xxpairs), 
								stringsAsFactors = FALSE
							)
	prog.bar.rdp.extm <- progress::progress_bar$new(total = nrow(xxpairs.conv.rev.raw))
	prog.bar.rdp.extm$tick(0)
	xxpairs.conv.rev.ext <- apply(xxpairs.conv.rev.raw, MARGIN = 1, function(x) {
		if (is.na(x[1])) {
			res <- c(x[2], x[1])
		} else {
			res <- if (x[1] > x[2]) c(x[1], x[2]) else c(x[2], x[1])
		}
		prog.bar.rdp.extm$tick()
		res
	})
	conv.rev.mat <- matrix(xxpairs.conv.rev.ext, ncol = 2, byrow = TRUE)
	xxpairs.rev.dup <- xxpairs[conv.rev.mat[,1], ]
	xxpairs.final <- unique(xxpairs.rev.dup)
	# return
	xxpairs.final
}


# bbpairs.final <- LREMapPairs(bbpairs, ligand.m.db, receptor.m.db, ecm.m.db)
# bbpairs.final <- RemoveDupPairs(bbpairs.final)
# done, with 45417 rows, at 2019.10.26

# dicide whether using anotation or not
ifUseAnnotation <- TRUE
#
#sspairs.final <- DoPartUnique(sspairs.recalc, c(1:2))
# make sure the type is right
sspairs.recalc[, 1] <- as.integer(sspairs.recalc[, 1])
sspairs.recalc[, 2] <- as.integer(sspairs.recalc[, 2])
#
sspairs.final <- FastAlignPairs(sspairs.recalc, 8)
sspairs.final <- sspairs.final[order(sspairs.final[, 3], sspairs.final[, 4],  # geneA geneB
	sspairs.final[, 12], decreasing = c(FALSE, FALSE, TRUE), method = "radix"), ]  # score as sub order
sspairs.final <- DoPartUnique(sspairs.final, 3:4)
#
# 5758309 rows at 2020.12.20 Human STRING v11
#
sspairs.final <- LREMapPairs(sspairs.final, ligand.anno.m.db, receptor.anno.m.db, ecm.anno.m.db, ifanno = ifUseAnnotation)
#sspairs.final <- RemoveDupPairs(sspairs.final)
# rearrange cols
if (ifUseAnnotation) {
	sspairs.runout <- sspairs.final[, c(1:8, 20:21, 14:17, 9:11)]
} else {
	sspairs.runout <- sspairs.final[, c(1:8, 14:15, 9:11)]
}
# check
str(sspairs.runout)
# done, with 42919 rows, at 2019.11.15, all mt-DNAs are excluede due to LRE list
# done, with 42919 rows, at 2019.11.09
#
# MORE TO SAY
# Npnt exists in both Ligand, Receptor, ECM
# while it is annotated to be in Extracelluar Space in Ligand, Receptor
# and be in Plasma Membrane in ECM
# ----
# 2019.11.09, concerning it may be marked wrong in ECM, which is contradictory
# As function uses Ligand as prior
# So, this warning is just ignored.
#

#aapairs.final <- DoPartUnique(aapairs.recalc, c(1:2))
aapairs.final <- FastAlignPairs(aapairs.recalc, 8)
aapairs.final <- aapairs.final[order(aapairs.final[, 3], aapairs.final[, 4],  # geneA geneB
	aapairs.final[, 12], decreasing = c(FALSE, FALSE, TRUE), method = "radix"), ]  # score as sub order
aapairs.final <- DoPartUnique(aapairs.final, 3:4)
aapairs.final <- LREMapPairs(aapairs.final, ligand.anno.m.db, receptor.anno.m.db, ecm.anno.m.db, ifanno = ifUseAnnotation)
#aapairs.final <- RemoveDupPairs(aapairs.final)
# rearrange cols
if (ifUseAnnotation) {
	aapairs.runout <- aapairs.final[, c(1:8, 20:21, 14:17, 9:11)]
} else {
	aapairs.runout <- aapairs.final[, c(1:8, 14:15, 9:11)]
}
# done, with 3311 rows, at 2019.10.26
#
#
#
# further screen
# remove <701 predicted score & > 700 known score in aapairs.runout, there are several.
# remove <701 known score in sspairs.runout, actually is 0
# !!!!!!  
# for these part, just done in Xlsx. Have not directly in codes.
sspairs.runout.ext
aapairs.runout.ext
#
# write.csv(ggpairs.final, file="LREpairs.csv", row.names=FALSE)








# [METHOD]
# @param ncols. integer. specify the ncol() of something
# e.g.
# if ncols = 4, returns c(2,1,4,3)
#
ReverseOddEvenCols <- function(ncols) {
	len.all <- ncols %/% 2
	val.even <- 2 * 1:len.all
	val.odds <- (2 * 1:len.all) -1
	new.serial <- NULL
	for (i in 1:len.all) {
		new.serial <- c(new.serial, val.even[i], val.odds[i])
	}
	if (ncols %% 2 != 0) {  # odds cols
		# leave last one unchanged
		new.serial <- c(new.serial, ncols)
	}
	# return
	new.serial
}
# [METHOD]
# to rebuild LREval in xxpairs
LRERecalcValPairs <- function(xxpairs) {
 LRE.bin.list <- c(
	NA, # 0
	"ECM", 	# 1
	"Receptor", # 2
	"Receptor, ECM", # 3
	"Ligand", # 4
	"Ligand, ECM", # 5
	"Ligand, Receptor", # 6
	"Ligand, Receptor, ECM" # 7
 )
 LREval.A <- match(xxpairs$inter.LREtype.A, LRE.bin.list)
 LREval.B <- match(xxpairs$inter.LREtype.B, LRE.bin.list)
 xxpairs$inter.LREval.A <- LREval.A
 xxpairs$inter.LREval.B	<- LREval.B
 # return
 xxpairs
}
# [METHOD] group
# for find l-l, l-r, l-e, ..., all kinds of interaction pairs
#
# 
# @param xxpairs format must be like bbpairs.final or aapairs.final
# @param inter.type1 character. only in "ligand", "receptor", "ecm"
# @param inter.type2 character. only in as inter.type1
# @param orderBywhich integer. if odds by inter.A, else by inter.B
#
# @import ReverseOddEvenCols
# 
LREFindDesiredPairs <- function(xxpairs, inter.type1, inter.type2, orderBywhich=1) {
	# ligand 	4,5,6,7
	LREval.ligand	<- c(4,5,6,7)
	# receptor 	2,3,6,7
	LREval.receptor <- c(2,3,6,7)
	# ecm 		1,3,5,7
	LREval.ecm 		<- c(1,3,5,7)
	# selection
	undefined.input.isornot <- FALSE
	val.sel1 <- LREval.receptor  # if select "receptor" or other undefined value
	if (inter.type1 == "ligand") val.sel1 <- LREval.ligand else {
		if (inter.type1 == "ecm") val.sel1 <- LREval.ecm else {
			if (inter.type1 != "receptor") undefined.input.isornot <- TRUE
		}
	}
	val.sel2 <- LREval.ligand    # if select "ligand" or other undefined value
	if (inter.type2 == "receptor") val.sel2 <- LREval.receptor else {
		if (inter.type2 == "ecm") val.sel2 <- LREval.ecm else {
			if (inter.type2 != "ligand") undefined.input.isornot <- TRUE
		}
	}
	if (undefined.input.isornot) { # if undefined input exists
		val.sel1 <- LREval.receptor
		val.sel2 <- LREval.ligand
		warning("Using default selection, return Ligand-Receptor pairs.")
	}
	# get conserved right order pairs
	as.conv.A <- which(xxpairs$inter.LREval.A %in% val.sel1)
	as.conv.B <- which(xxpairs$inter.LREval.B %in% val.sel2)
	as.pairs.conv <- intersect(as.conv.A, as.conv.B)
	# get reversed right order pairs
	as.rev.B <- which(xxpairs$inter.LREval.B %in% val.sel1)
	as.rev.A <- which(xxpairs$inter.LREval.A %in% val.sel2)
	as.pairs.rev <- intersect(as.rev.B, as.rev.A)
	# get pairs of conv , rev
	rpairs.conv <- xxpairs[as.pairs.conv, ]
	rpairs.rev  <- xxpairs[as.pairs.rev, ReverseOddEvenCols(ncol(xxpairs))]
	colnames(rpairs.rev) <- colnames(xxpairs)
	# rbind
	rpairs.result <- rbind(rpairs.conv, rpairs.rev)
	rpairs.result <- unique(rpairs.result)
	if (orderBywhich %% 2 != 0) {
		rpairs.result <- rpairs.result[order(rpairs.result$inter.GeneName.A), ]
	} else {
		rpairs.result <- rpairs.result[order(rpairs.result$inter.GeneName.B), ]
	}
	# return
	rpairs.result
}




