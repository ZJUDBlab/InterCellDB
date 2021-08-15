# test database overlaps and confidence
# 2020.11.28 ver1

library(CellTalkDB)  # load CellTalkDB database
library(nichenetr)  # load NicheNet database
library(dplyr)

cross.GetStDB <- function(
	user.db,
	ref.taxonomy = "human",
	colnames.used.from = "from",
	colnames.used.to = "to",
	...
) {
	if (ref.taxonomy == "human") {
		this.ref.db <- genes.human.ref.db
	} else {
		if (ref.taxonomy == "mouse") {
			this.ref.db <- genes.mouse.ref.db
		} else {
			stop("No ref database!")
		}
	}

	inter.from <- user.db[, colnames.used.from, drop = FALSE]
	inter.from <- cbind(inter.from, row.id = 1:nrow(inter.from))
	inter.to <- user.db[, colnames.used.to, drop = FALSE]
	inter.to <- cbind(inter.to, row.id = 1:nrow(inter.to))

	new.inter.from <- Tool.AddUserRestrictDB(inter.from, this.ref.db)
	new.inter.to <- Tool.AddUserRestrictDB(inter.to, this.ref.db)
	new.inter.from <- new.inter.from[order(new.inter.from$row.id), ]
	new.inter.to <- new.inter.to[order(new.inter.to$row.id), ]
	colnames(new.inter.from)[1:2] <- c("inter.GeneID.A", "inter.GeneName.A")
	colnames(new.inter.to)[1:2] <- c("inter.GeneID.B", "inter.GeneName.B")

	new.inter.db <- left_join(new.inter.from, new.inter.to)
	new.inter.db$row.id <- NULL
	new.inter.db <- cbind(new.inter.db, user.db[, setdiff(colnames(user.db), c(colnames.used.from, colnames.used.to))])
	if (ncol(new.inter.db) == 4) {
		new.inter.db <- new.inter.db[, c(1,3,2,4)]  # align like CellTalkDB
	} else {
		new.inter.db <- new.inter.db[, c(1,3,2,4, 5:ncol(new.inter.db))]
	}

	return(new.inter.db)
}

cross.CheckOverlap <- function(
	check.target.db,
	check.ref.db,
	target.name = "NicheNet",
	ref.name = "CellTalkDB",
	check.cols = 4,
	...
) {
	to.check.target <- cbind(check.target.db, "cross.xxx.source" = target.name, "cross.xxx.row.id" = 1:nrow(check.target.db))
	to.check.ref <- cbind(check.ref.db, "cross.xxx.source" = ref.name, "cross.xxx.row.id" = 1:nrow(check.ref.db))
	check.merge.db <- rbind(to.check.ref[, c(1:check.cols, (ncol(to.check.ref) - 1):ncol(to.check.ref))],
		to.check.target[, c(1:check.cols, (ncol(to.check.target) - 1):ncol(to.check.target))])
	check.merge.db <- DoPartUnique(check.merge.db, c(check.cols - 1, check.cols))
	#
	check.rest.ids <- check.merge.db[which(check.merge.db[, "cross.xxx.source"] == target.name), "cross.xxx.row.id"]
	overlap.cnt <- nrow(to.check.target) - length(check.rest.ids)
	return(list(cnt = overlap.cnt, rest.db = to.check.target[check.rest.ids, 1:(ncol(to.check.target) - 2)]))
}

# process
# readin raw db
n2015.lrdb <- readRDS("../2015Nature-LR-database.rds")
n2020.lrdb <- readRDS("../2020Nature-upd-LR-database.rds")
cp.lrdb <- readRDS("../CellPhoneDB-gene-pairs.rds")

##
## [TO NOTE]: LR db has the direction of Ligand->Receptor. Here, I remove this direction.
## The process followed only compares the unique gene pairs that database presents.
##

# 2015 Nature LR db
new.db.2015 <- cross.GetStDB(n2015.lrdb, "human", "Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")
new.db.2015 <- DoPartUnique(FastAlignPairs(new.db.2015, 4), 3:4)  # 2557-> 2555 rows

# 2020 Nature LR db
new.db.2020 <- cross.GetStDB(n2020.lrdb, "human", "Ligand.gene.symbol", "Receptor.gene.symbol")
new.db.2020 <- DoPartUnique(FastAlignPairs(new.db.2020, 4), 3:4)  # 2293-> 2248 rows

# CellPhoneDB transed LR db
new.cpdb <- cross.GetStDB(cp.lrdb, "human", "genename.A", "genename.B")
new.cpdb <- DoPartUnique(FastAlignPairs(new.cpdb, 4), 3:4)  # 2789-> 1580 rows

# NicheNet lr_network
new.db.nichenet <- cross.GetStDB(lr_network, "human", "from", "to")
new.db.nichenet <- FastAlignPairs(new.db.nichenet, 4)
# levels(factor(new.db.nichenet$database))
tmp.order.sub.db <- c("ramilowski", "kegg", "guide2pharmacology", "ppi_prediction", "ppi_prediction_go")
tmp.sub.list <- integer()
for (i in tmp.order.sub.db) {
	tmp.sub.list <- c(tmp.sub.list, which(new.db.nichenet$database == i))
}
new.db.nichenet <- new.db.nichenet[tmp.sub.list, ]
new.db.nichenet <- DoPartUnique(new.db.nichenet, 3:4) 
# total 11962 rows
# get ramilowski 1164, guide2pharmacology 24, kegg 191, ppi_* takes the rest 6280 and 4303.

# String all DB
new.db.string <- readRDS("../res-db-packed/STRING-9606-20201204.rds")
new.db.string <- new.db.string[, c(1:4, 9:12)]  # may check
#new.db.string <- cross.GetStDB(string.all.db, "human", "inter.GeneName.A", "inter.GeneName.B")



##### other part by the way
# test data loss - 2015 Nature LR db
new.db.2015.t1 <- cross.GetStDB(n2015.lrdb, "human", "Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")
new.db.2015.t1 <- FastAlignPairs(new.db.2015.t1, 4)
#
new.db.2015.t2 <- DoPartUnique(new.db.2015.t1, 3:4)
overlap.2015.rownames <- setdiff(rownames(new.db.2015.t1), rownames(new.db.2015.t2))
overlap.2015.ids <- which(rownames(new.db.2015.t1) %in% overlap.2015.rownames)
# finally catch the duplicates in final database,
# 1. AREG~EGFR 2. AREG~ERBB3
# the corresponding ones in original database are listed as follows,
# 1. AREG~EGFR 2. AREG~ERBB3 3. AREGB~EGFR 4. AREGB~ERBB3



###### main process in cross-comparison
### [TO NOTE] as all data has been FastAlign, *.rev will be 0, to be as the checkpoint
### part 1: use 2015 Nature as ref
ref.rami.conv <- paste(new.db.2015[, "inter.GeneName.A"], new.db.2015[, "inter.GeneName.B"], sep = "~")
ref.rami.rev  <- paste(new.db.2015[, "inter.GeneName.B"], new.db.2015[, "inter.GeneName.A"], sep = "~")
# [TO NOTE] intersect the above 2 (conv & rev) gives NO overlaps.

### part 1.ex: use 2020 Nature as ref
ref.rami.conv <- paste(new.db.2020[, "inter.GeneName.A"], new.db.2020[, "inter.GeneName.B"], sep = "~")
ref.rami.rev  <- paste(new.db.2020[, "inter.GeneName.B"], new.db.2020[, "inter.GeneName.A"], sep = "~")
# [TO NOTE] intersect the above 2 (conv & rev) gives 21 overlaps.
# BUT, these are all A-A the gene pairs involving only 1 gene itself.

#
## part 1.1: take NicheNet and CellTalkDB in
# NicheNet lr_network
tg1.1.niche <- paste(new.db.nichenet[, "inter.GeneName.A"], new.db.nichenet[, "inter.GeneName.B"], sep = "~")
inds.tg1.1.niche.conv <- which(tg1.1.niche %in% ref.rami.conv)
inds.tg1.1.niche.rev <- which(tg1.1.niche %in% ref.rami.rev)
if (length(inds.tg1.1.niche.rev) != 0) { stop("Unexpected niche x rami.rev is not 0.") }
inds.tg1.1.niche <- c(inds.tg1.1.niche.conv, inds.tg1.1.niche.rev)
tg1.1.niche.x.ref <- unique(tg1.1.niche[inds.tg1.1.niche])
# total 1656 matches in 2015 Nature
# total 1513 matches in 2020 Nature

## not main process
# check matches in NicheNet origins
tg1.1.niche.matches.group <- tapply(1:length(inds.tg1.1.niche), new.db.nichenet[inds.tg1.1.niche, "database"], length)
#@# matches that database noted as ramilowski: 1163, compare to total ramilowski 2015, percentage: 0.455
# - for the 1 not matched, ramilowski in lr_network is 1164 records
inds.tg1.1.niche.is.rami <- intersect(which(new.db.nichenet$database == "ramilowski"), inds.tg1.1.niche)
inds.tg1.1.niche.not.match.in.rami <- setdiff(which(new.db.nichenet$database == "ramilowski"), inds.tg1.1.niche.is.rami)
# - finally catch the crime. CCL3L3~CCR1 that is not recorded by ramilowski in 2015 lrdb.
#@# matches ignoring the manually recorded as ramilowski: 1656, compare to total ramilowski 2015, percentage: 0.648


# CellTalkDB
tg1.1.celltalk <- paste(new.db.string[, "inter.GeneName.A"], new.db.string[, "inter.GeneName.B"], sep = "~")
inds.tg1.1.celltalk.conv <- which(tg1.1.celltalk %in% ref.rami.conv)
inds.tg1.1.celltalk.rev <- which(tg1.1.celltalk %in% ref.rami.rev)
if (length(inds.tg1.1.celltalk.rev) != 0) { stop("Unexpected celltalk x rami.rev is not 0.") }
inds.tg1.1.celltalk <- c(inds.tg1.1.celltalk.conv, inds.tg1.1.celltalk.rev)
tg1.1.celltalk.x.ref <- unique(tg1.1.celltalk[inds.tg1.1.celltalk])
# total 2442 matches in 2015 Nature
# total 2140 matches in 2020 Nature

# intersection NicheNet and CellTalkDB
tg1.1.niche.rev <- paste(new.db.nichenet[, "inter.GeneName.B"], new.db.nichenet[, "inter.GeneName.A"], sep = "~")
inds.tg1.1.celltalk.x.niche.conv <- which(tg1.1.celltalk %in% tg1.1.niche)
inds.tg1.1.celltalk.x.niche.rev <- which(tg1.1.celltalk %in% tg1.1.niche.rev)
if (length(inds.tg1.1.celltalk.x.niche.rev) != 0) { stop("Unexpected celltalk x niche.rev is not 0.") }
inds.tg1.1.celltalk.x.niche <- c(inds.tg1.1.celltalk.x.niche.conv, inds.tg1.1.celltalk.x.niche.rev)
tg1.1.celltalk.x.niche <- unique(tg1.1.celltalk[inds.tg1.1.celltalk.x.niche])
# total 11256 matches for celltalk x niche

# intersection of 3 DBs, NicheNet, CellTalkDB, 2015/2020 Nature.
tg1.1.3db.overlap <- intersect(ref.rami.conv, intersect(tg1.1.niche, tg1.1.celltalk))
# total 1624 matches for 2015 Nature 3 overlap
# total 1488 matches for 2020 Nature 3 overlap

## part 1.2: take CellPhoneDB and CellTalkDB in
# CellPhoneDB
tg1.2.cpdb <- paste(new.cpdb[, "inter.GeneName.A"], new.cpdb[, "inter.GeneName.B"], sep = "~")
inds.tg1.2.cpdb.conv <- which(tg1.2.cpdb %in% ref.rami.conv)
inds.tg1.2.cpdb.rev <- which(tg1.2.cpdb %in% ref.rami.rev)
if (length(inds.tg1.2.cpdb.rev) != 0) { stop("Unexpected cpdb x rami.rev is not 0.") }
inds.tg1.2.cpdb <- c(inds.tg1.2.cpdb.conv, inds.tg1.2.cpdb.rev)
tg1.2.cpdb.x.ref <- unique(tg1.2.cpdb[inds.tg1.2.cpdb])
# total 936 matches in 2015 Nature
# total 980 matches in 2020 Nature

# CellTalkDB, the result same as the above
# tg1.1.celltalk, tg1.1.celltalk.x.ref
tg1.2.celltalk <- tg1.1.celltalk
tg1.2.celltalk.x.ref <- tg1.1.celltalk.x.ref

# intersection CellPhoneDB and CellTalkDB
tg1.2.cpdb.rev <- paste(new.cpdb[, "inter.GeneName.B"], new.cpdb[, "inter.GeneName.A"], sep = "~")
inds.tg1.2.celltalk.x.cpdb.conv <- which(tg1.2.celltalk %in% tg1.2.cpdb)
inds.tg1.2.celltalk.x.cpdb.rev <- which(tg1.2.celltalk %in% tg1.2.cpdb.rev)
if (length(inds.tg1.2.celltalk.x.cpdb.rev) != 0) { stop("Unexpected celltalk x cpdb.rev is not 0.") }
inds.tg1.2.celltalk.x.cpdb <- c(inds.tg1.2.celltalk.x.cpdb.conv, inds.tg1.2.celltalk.x.cpdb.rev)
tg1.2.celltalk.x.cpdb <- unique(tg1.2.celltalk[inds.tg1.2.celltalk.x.cpdb])
# total 1401 matches in celltalk x cpdb

# intersection of 3 DBs, CellPhoneDB, CellTalkDB, 2015/2020 Nature.
tg1.2.3db.overlap <- intersect(ref.rami.conv, intersect(tg1.2.cpdb, tg1.2.celltalk))
# total 929 matches for 2015 Nature 3 overlap
# total 968 matches for 2020 Nature 3 overlap

