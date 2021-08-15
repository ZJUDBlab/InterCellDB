
# -------------------------------------
# check [MANUAL] tags in this file, 
# check [TODO] tags as well.

# -------------------------------------


#### reference database
# entrez database
entrez.9606.db  # use standard merge database at 2020.03.22
entrez.10090.db  # use standard merge database at 2020.03.22

# keyword list
uniprot.keywords.all <- read.table("../keywords-all.tab.txt", 
		header = TRUE,
		sep = "\t",
		quote = "",
		comment.char = "",
		stringsAsFactors = FALSE)
# remove the column "Description", and saveRDS



#### process database
# 9606 uniprot genes - all
uniprot.9606.genes.all <- read.table("../DB-based-new/uniprot-human-20200322.txt",
		header = TRUE,
		sep = "\t",
		quote = "",
		comment.char = "",
		stringsAsFactors = FALSE)
# second time, use RDS "../DB-based-new/uniprot-9606-genes-all-20200227.rds"

# 10090 uniprot genes -all
uniprot.10090.genes.all <- read.table("../DB-based-new/uniprot-mouse-20200322.txt",
		header = TRUE,
		sep = "\t",
		quote = "",
		comment.char = "",
		stringsAsFactors = FALSE)

# ------

# entrez.9606.db <- readRDS("../DB-based-new/entrez-9606-20200226.rds")
uniprot.9606.genes.all <- readRDS("../DB-based-new/uniprot-9606-genes-all-20200227.rds")
uniprot.keywords.simple <- readRDS("../uniprot-keywords-simple-20200226.rds")

# [MANUAL] [TO change species] get all below variables reset
# -------------------------------------
entrez.ref.db <- entrez.9606.db
uniprot.tax.genes.all <- uniprot.9606.genes.all
uniprot.keywords.ref <- uniprot.keywords.simple

colnames(uniprot.tax.genes.all)[10] <- "Cross.reference.GeneID"

# -------------------------------------



### map every entry to NCBI gene database
map.subset.genes.all <- uniprot.tax.genes.all[, c("Entry", "Gene.names", "Cross.reference.GeneID")]
# save subset NA
map.subset.genes.cross.NAs <- map.subset.genes.all[which(map.subset.genes.all$Cross.reference.GeneID == ""), ]
# remove NAs
map.subset.genes.all <- map.subset.genes.all[-which(map.subset.genes.all$Cross.reference.GeneID == ""), ]  # 18858 rows at 2020.02.27

## get those with only 1 or multiply mapped GeneID apart
map.subset.genes.all$Cross.reference.GeneID <- sapply(map.subset.genes.all$Cross.reference.GeneID, function(x) {
	y <- substring(x, 1, nchar(x) - 1)
	as.character(y)
	})
# check ';' to figure out those with multiply mapped GeneIDs
inds.with.mul.map <- grep(";", map.subset.genes.all[, "Cross.reference.GeneID"], fixed = TRUE)
#map.subset.genes.onemap <- map.subset.genes.all
#map.subset.genes.multimap <- NULL
map.subset.genes.onemap <- map.subset.genes.all[-inds.with.mul.map, ]
map.subset.genes.multimap <- map.subset.genes.all[inds.with.mul.map, ]


# process for determination of the actual mapping genes in *.multimap subset
if (!is.null(map.subset.genes.multimap) && nrow(map.subset.genes.multimap) > 0) {
	## 1st step, split the multi gene.names
	inds.multi.gn <- grep(";", map.subset.genes.multimap[, "Gene.names"], fixed = TRUE)
	to.repack.single.gn <- map.subset.genes.multimap[-inds.multi.gn, ]
	to.repack.multi.gn <- map.subset.genes.multimap[inds.multi.gn, ]
	# 1.1 step, process the single.gene.names
	if (nrow(to.repack.single.gn) > 0) {
		tmp.new.single.gn <- apply(to.repack.single.gn, MARGIN = 1, FUN = function(x) {
			tmp.m.IDs <- strsplit(x["Cross.reference.GeneID"], split = ";", fixed = TRUE)[[1]]
			data.frame(Entry = as.character(x["Entry"]), Gene.names = as.character(x["Gene.names"]),
				Cross.reference.GeneID = tmp.m.IDs, stringsAsFactors = FALSE)
		})
		tmp.new.df <- data.frame()
		for (tmp.i in tmp.new.single.gn) {
			tmp.new.df <- rbind(tmp.new.df, tmp.i)
		}
		# set back
		tmp.save.colnames <- colnames(to.repack.single.gn)
		to.repack.single.gn <- tmp.new.df
		colnames(to.repack.single.gn) <- tmp.save.colnames
	}
	# 1.2 step, process the multi.gene.names
	if (nrow(to.repack.multi.gn) > 0) {
		tmp.new.multi.gn <- apply(to.repack.multi.gn, MARGIN = 1, FUN = function(x) {
			tmp.m.gns <- strsplit(x["Gene.names"], split = "; ", fixed = TRUE)[[1]]
			tmp.m.IDs <- strsplit(x["Cross.reference.GeneID"], split = ";", fixed = TRUE)[[1]]
			if (length(tmp.m.gns) != length(tmp.m.IDs)) {
				#warning(paste0("The length of 'multi.gene.names' is NOT EQUAL to 'Cross.reference.GeneID'!\n", 
				#	as.character(x["Entry"])))
				# note: never mind if not equal, as using different strategy to treat them, see below.
			}
			# Uniprot don't make the gene.names one-by-one mapping to the following GeneIDs
			# So, pack up all gene names
			tmp.list.gns <- sapply(tmp.m.gns, function(x) {
				tmp.split <- strsplit(x, split = " ", fixed = TRUE)[[1]]
				as.character(tmp.split)
			})
			tmp.vec.gns <- as.character(unlist(tmp.list.gns))
			data.frame(Entry = as.character(x["Entry"]), 
				Gene.names = paste0(tmp.vec.gns, collapse = " "), 
				Cross.reference.GeneID = tmp.m.IDs, stringsAsFactors = FALSE)
		})
		tmp.new.df <- data.frame()
		for (tmp.i in tmp.new.multi.gn) {
			tmp.new.df <- rbind(tmp.new.df, tmp.i)
		}
		# set back
		tmp.save.colnames <- colnames(to.repack.multi.gn)
		to.repack.multi.gn <- tmp.new.df
		colnames(to.repack.multi.gn) <- tmp.save.colnames
	}
	# 2nd done step, merge single & multi, and set back
	map.subset.genes.multimap <- rbind(to.repack.single.gn, to.repack.multi.gn)
}

# merge all
adj.subset.genes.all <- rbind(map.subset.genes.onemap, map.subset.genes.multimap, map.subset.genes.cross.NAs)
	#*!*# for automatic process, recheck and set stop point
	inds.auto.check.onegene <- grep(";", adj.subset.genes.all[, "Cross.reference.GeneID"], fixed = TRUE)
	if (length(inds.auto.check.onegene) > 0) {
		stop("Auto remove multi-gene-mapping failed! Please recheck this process!")
	}



# map back to NCBI gene database
prog.bar.adj.e.m.t <- progress::progress_bar$new(total = nrow(adj.subset.genes.all))
prog.bar.adj.e.m.t$tick(0)
adj.entrez.map.tres <- apply(adj.subset.genes.all, MARGIN = 1, entrez.ref.db = entrez.ref.db,
	function(x, entrez.ref.db) {
		tmp.GeneID <- as.integer(x["Cross.reference.GeneID"])
		tref.genenames.list <- strsplit(as.character(x["Gene.names"]), split = " ", fixed = TRUE)[[1]]
		# do map to NCBI gene based on GeneID
		ind.match <- which(entrez.ref.db[, "GeneID"] == tmp.GeneID)
		if (length(ind.match) >= 2) {
			stop("Unexpected duplicate matches. Please recheck database if it was broken by manual operations!")
		}
		t.res.match.GeneID <- NA
		t.res.checkname <- NA  # if checkname is TRUE, "CHECK.TRUE" is set, else the 
		# 1st, Authorized gene name
		if (length(ind.match) > 0) {
			tm.auth.genename <- entrez.ref.db[ind.match, "Symbol_from_nomenclature_authority"]
			if (tm.auth.genename %in% tref.genenames.list) {
				t.res.match.GeneID <- tmp.GeneID
				t.res.checkname <- "CHECK.TRUE"
			} else {  # use symbol and synonyms
				tm.other.genename <- c(entrez.ref.db[ind.match, "Symbol"], strsplit(entrez.ref.db[ind.match, "Synonyms"], split = "|", fixed = TRUE)[[1]])
				if (sum(tref.genenames.list %in% tm.other.genename) > 0) {
					t.res.match.GeneID <- tmp.GeneID
					t.res.checkname <- "SYN.MAP.TRUE"
				}
			}
			# else leave out NAs
		}
		# 2nd, try to use uniprot.gene.names to get matches
		if (is.na(t.res.match.GeneID)) {
			ind.name.tries <- integer()
			for (i.name in tref.genenames.list) {
				ind.tname.m <- which(entrez.ref.db[, "Symbol_from_nomenclature_authority"] == i.name)
				ind.name.tries <- c(ind.name.tries, ind.tname.m)
			}
			ind.name.tries <- unique(ind.name.tries)
			if (length(ind.name.tries) > 1) {
				warning(paste0("Uniprot gene.names map to different NCBI gene: ", paste0(tref.genenames.list, collapse = " ")),
					paste0(" Overwrite with the 1st match: ", tref.genenames.list[1]))
				ind.name.tries <- ind.name.tries[1]  # be overwritten here
			}
			if (length(ind.name.tries) == 1) {
				t.res.match.GeneID <- entrez.ref.db[ind.name.tries, "GeneID"]
				t.res.checkname <- "REMAP.AUTH.TRUE"
			}
		}
		prog.bar.adj.e.m.t$tick()
		c(as.character(t.res.match.GeneID), as.character(t.res.checkname))
	}
)
adj.entrez.map.tres <- data.frame(t(adj.entrez.map.tres), stringsAsFactors = FALSE)
colnames(adj.entrez.map.tres) <- c("map.GeneID", "check.res")
adj.entrez.map.tres[, 1] <- as.integer(adj.entrez.map.tres[, 1])
#
# ! Use a small part of manual check on unmatched ones
# ! Most unmatched rows contain informations out-of-date or not included in this release of entrez.ref.db.
# [TODO] a more elegant way to deal with it
#
ind.unmatched.NAs <- which(is.na(adj.entrez.map.tres[, "map.GeneID"]))
leavout.items <- adj.subset.genes.all[ind.unmatched.NAs, ]  # leave out those with no explict GeneID
# get not NA rows
adj.genes.result.b1 <- cbind(adj.subset.genes.all, adj.entrez.map.tres)
if (length(ind.unmatched.NAs) > 0) {
	adj.genes.result.b1 <- adj.genes.result.b1[-ind.unmatched.NAs, ]
} else {
	adj.genes.result.b1 <- adj.genes.result.b1
}
# get result: adj.genes.result.b1
nrow(adj.genes.result.b1)  # 20455 rows at 2020.02.29
# 20456 rows at 2020.12.20, for human
# 16832 rows at 2020.12.20, for mouse


if (FALSE) {  # deprecated at 2020.03.12, no meaning to check duplicate entries, only need to map it to EntrezID and do unique
	# check one more time, as multi* split genenames & cross.ref.s to make duplicate entries
	adj.genes.recheck <- adj.genes.result.b1
	avec.recheck.entry <- factor(adj.genes.recheck[, 1])
	avec.find.dup <- tapply(rep(1, length(avec.recheck.entry)), avec.recheck.entry, length)
	avec.dup.entries <- names(avec.find.dup)[which(avec.find.dup > 1)]
	#
	# recorded at 2020.02.29
	# how to distinguish these duplicate entries
	# 1. Complex. To recognize this, current strategy is 
	#		all names(2:<last>) can be grep by the first one by ^<first>
	#		and if not complex,
	#			do further manual check
	# 2. Single recheck. NOTE: Real duplicate or ERROR when only single gene name appears. 
	#    To deal with this, current strategy is
	#		auto-check among each other, 
	#			0. all duplicate ones are the same, save only one of those, and leave others out.
	#			and if any conflicts exist, 
	#			1. save once, Use the gene name to match to NCBI gene database again, and get the matched GeneID as final result
	#			2. if step-1 failed, do some manual works to check through
	#			3. all above steps are failed, just leave it out
	# ex. Unknown. To recognize this, do some manual works 
	#
	# !NOTE! [MANUAL]
	# 	at 2020.02.29 version of result, all manual works above will be seen as leave-outs.
	# -----------------------------
	# 
	adj.genes.result.c2 <- adj.genes.recheck
	adj.c2.feedback <- list()
	# process for these duplicate entries
	if (length(avec.dup.entries) == 0) {
		# no entry is duplicate, GOOD!
	} else {  # allow modifcation direcly to adj.genes.result.c2
		# get all split entries of Complex through, and leave out the rest unknown ones
		ind.dup.matchback <- integer(0)
		ind.dup.leavout <- integer(0)
		reason.dup.leaveout <- character(0)  # one-by-one matched with the ind.dup.leavout
		for (i.dup in avec.dup.entries) {
			ind.this.dup <- which(adj.genes.recheck[, "Entry"] == i.dup)  # implicitly > 1
			# get gene.names
			this.genenames <- strsplit(adj.genes.recheck[ind.this.dup[1], "Gene.names"], split = " ", fixed = TRUE)[[1]]
			this.ref.name <- this.genenames[1]
			this.other.name <- this.genenames[-1]
			if (length(this.other.name) > 0) {
				# judge Complex process
				ind.to.ref <- grep(paste0("^", this.ref.name), this.other.name)
				if (length(ind.to.ref) == length(this.other.name)) {
					ind.dup.matchback <- c(ind.dup.matchback, ind.this.dup)  # get all through
				} else {
					ind.dup.leavout <- c(ind.dup.leavout, ind.this.dup)
					reason.dup.leaveout <- c(reason.dup.leaveout, rep("NOT.COMPLEX", times = length(ind.this.dup)))
				}
			} else {
				# check if conflicts exist on single gene name
				this.map.GeneIDs <- adj.genes.recheck[ind.this.dup, "map.GeneID"]
				this.map.GeneIDs <- unique(this.map.GeneIDs)
				if (length(this.map.GeneIDs) == 1) {  # strategy: 2-step 0
					ind.dup.matchback <- c(ind.dup.matchback, ind.this.dup[1])  # remove all other duplicate ones, because they are the same
					ind.dup.leavout <- c(ind.dup.leavout, ind.this.dup[-1])
					reason.dup.leaveout <- c(reason.dup.leaveout, rep("REMOVE.DUP", times = length(ind.this.dup[-1])))
				} else {  # conflicts exist
					this.single.gene.name <- adj.genes.recheck[ind.this.dup[1], "Gene.names"]
					ind.save.match <- which(entrez.ref.db[, "Symbol_from_nomenclature_authority"] == this.single.gene.name)
					if (length(ind.save.match) > 0) {  # strategy: 2-step 1
						# direct modification
						adj.genes.result.c2[ind.this.dup[1], "map.GeneID"] <- entrez.ref.db[ind.save.match, "GeneID"]
						adj.genes.result.c2[ind.this.dup[1], "check.res"] <- "DUP.OVREWRITE"
						#
						ind.dup.matchback <- c(ind.dup.matchback, ind.this.dup[1])
						ind.dup.leavout <- c(ind.dup.leavout, ind.this.dup[-1])
						reason.dup.leaveout <- c(reason.dup.leaveout, rep("SAVE.LEAVES", times = length(ind.this.dup[-1])))
					} else {  # strategy: 2-step >2
						ind.dup.leavout <- c(ind.dup.leavout, ind.this.dup)
						reason.dup.leaveout <- c(reason.dup.leaveout, rep("NEED.MANUAL", times = length(ind.this.dup)))
					}
				}
			}
		}
		# contruct the result
		adj.genes.result.c2 <- adj.genes.result.c2[-ind.dup.leavout, ]  # [TODO] how to construct
		# feed back infos
		adj.c2.feedback <- list(matchback = ind.dup.matchback,
			leaveout = ind.dup.leavout, reason.leaveout = reason.dup.leaveout)
	}
	# get result: adj.genes.result.c2
	nrow(adj.genes.result.c2)  # 19994 rows at 2020.02.29
	#
}



adj.entry.remap.c2 <- adj.genes.result.b1[, c("Entry", "map.GeneID")]
#
#
#
# -------------------------------------

# keyword subset  "Molecular function"
uniprot.subset.mfunc <- uniprot.keywords.ref[which(uniprot.keywords.ref$Category == "Molecular function"), ]
# keyword subset "Cellular component"
uniprot.subset.ccomp <- uniprot.keywords.ref[which(uniprot.keywords.ref$Category == "Cellular component"), ]

# start from uniprot.tax.genes.all to split keywords
skey.subset.genes.all <- uniprot.tax.genes.all[, c("Entry", "Keyword.ID", "Keywords")]
prog.bar.skey.sp.l <- progress::progress_bar$new(total = nrow(skey.subset.genes.all))
prog.bar.skey.sp.l$tick(0)
skey.splits.list <- apply(skey.subset.genes.all, MARGIN = 1, 
	mfunc.ref = uniprot.subset.mfunc, ccomp.ref = uniprot.subset.ccomp,
	function(x, mfunc.ref, ccomp.ref) {
		# keyword split
		this.entry <- as.character(x["Entry"])
		this.keys <- strsplit(as.character(x["Keyword.ID"]), split = "; ", fixed = TRUE)[[1]]
		# mfunc process
		inds.mfunc <- match(this.keys, mfunc.ref[, "Keyword.ID"])
		this.ret.mfunc <- mfunc.ref[inds.mfunc[which(!is.na(inds.mfunc))], ]
		if (nrow(this.ret.mfunc) > 0) {
		this.ret.mfunc <- cbind(Entry = this.entry, this.ret.mfunc, stringsAsFactors = FALSE)
		} else {
			this.ret.mfunc <- NULL
		}
		# ccomp process
		inds.ccomp <- match(this.keys, ccomp.ref[, "Keyword.ID"])
		this.ret.ccomp <- ccomp.ref[inds.ccomp[which(!is.na(inds.ccomp))], ]
		if (nrow(this.ret.ccomp) > 0) {
			this.ret.ccomp <- cbind(Entry = this.entry, this.ret.ccomp, stringsAsFactors = FALSE)
		} else {
			this.ret.ccomp <- NULL
		}
		prog.bar.skey.sp.l$tick()
		list(mfunc = this.ret.mfunc, ccomp = this.ret.ccomp)
	}
)
# re-collect splits
library(dplyr)

skey.splits.plain.1.list <- unlist(skey.splits.list, recursive = FALSE, use.names = FALSE)
skey.splits.collect <- bind_rows(skey.splits.plain.1.list)
skey.mfunc.collect <- skey.splits.collect[which(skey.splits.collect$Category == "Molecular function"), ]
skey.ccomp.collect <- skey.splits.collect[which(skey.splits.collect$Category == "Cellular component"), ]

#
#
# -------------------------------------

# merge databases
col.names.entrez.trunc <- c("GeneID", "Symbol_from_nomenclature_authority")
mdat.mfunc.tmp <- left_join(adj.entry.remap.c2, skey.mfunc.collect)
mdat.mfunc.final <- left_join(mdat.mfunc.tmp, entrez.ref.db[, col.names.entrez.trunc], by = c("map.GeneID" = "GeneID"))


#mdat.ccomp.tmp <- left_join(adj.entry.remap.c2, skey.ccomp.collect)
#mdat.ccomp.final <- left_join(mdat.ccomp.tmp, entrez.ref.db[, col.names.entrez.trunc], by = c("map.GeneID" = "GeneID"))

# [TODO] consider re-ordering the NAs and other available contents

# added at 2020.03.12
mdat.mfunc.final.ext <- DoPartUnique(mdat.mfunc.final, c(1,2,3))
mdat.mfunc.final.ext <- mdat.mfunc.final.ext[, c(1:2, 6, 3:5)]



# added at 2020.03.22
Tc.Cap.simple <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
}

Tc.Cap.simple.vec <- function(to.cap.vec) {
	unlist(lapply(to.cap.vec, FUN = Tc.Cap.simple))
}
#
mdat.mfunc.final.ext <- mdat.mfunc.final.ext[, c(2,3,1,4:6)]
colnames(mdat.mfunc.final.ext) <- c("GeneID", "Gene.name", "Entry.Uniprot", "Keyword.ID", "Keyword.Name", "Category")
mdat.mfunc.final.ext[which(is.na(mdat.mfunc.final.ext$Keyword.Name)), "Keyword.Name"] <- "Other"
mdat.mfunc.final.ext$Keyword.Name <- Tc.Cap.simple.vec(mdat.mfunc.final.ext$Keyword.Name)




