
library(dplyr)

# how to use the strategy described in article:COMPARTMENTS
# use package "GO.db"
library(GO.db)

#
# readin NCBI gene
# 9606
entrez.9606.db <- readRDS("../../homo-9606-workflow/DB-based/homo-9606-entrez-20191114.rds")
# 10090
entrez.10090.db <- readRDS("../../mus-10090-workflow/DB-based/mus-10090-entrez-20191020.rds")

#
# readin ensembl
# 9606
ensembl.9606.db <- read.csv("../../homo-9606-workflow/STRING/ensembl-9606-20191115.csv", stringsAsFactors=FALSE)
ensembl.9606.db <- ensembl.9606.db[which(ensembl.9606.db$Protein.stable.ID != ""), ]
#ensembl.9606.db <- ensembl.9606.db[which(!is.na(ensembl.9606.db$NCBI.gene.ID)), ]
# 10090
ensembl.10090.db <- read.csv("../../mus-10090-workflow/STRING/10090-csv-mgi-mart_export.txt", stringsAsFactors=FALSE)
ensembl.10090.db <- ensembl.10090.db[which(ensembl.10090.db$Protein.stable.ID != ""), ]
#


# new strategy to use ensembl
ensembl.tax.db[, 4] <- NULL  # new way to remove it, for human NCBI.gene.ID, for mouse MGI.ID
ensembl.tax.db <- unique(ensembl.tax.db)  # 110788 rows in *9606
tmp.check.cnt.p <- tapply(rep(1, times = nrow(ensembl.tax.db)), ensembl.tax.db$Protein.stable.ID, sum)
tmp.check.cnt.g <- tapply(rep(1, times = nrow(ensembl.tax.db)), ensembl.tax.db$Gene.stable.ID, sum)
tmp.check.cnt.gn <- tapply(rep(1, times = nrow(ensembl.tax.db)), ensembl.tax.db$Gene.name, sum)
tmp.len.dup <- length(which(tmp.check.cnt.p > 1)) + length(which(tmp.check.cnt.g > 1)) + length(which(tmp.check.cnt.gn > 1))
if (length(which(tmp.check.cnt.p > 1)) > 0) {
	stop("Proteins are duplicate!")
}
#

# compartment part
compt.template.colnames <- c("Protein.stable.ID", "Gene.name", "GO.ID", "GO.term", 
							 "Source", "Evidence", "score")
# -db_1- human knowledge location database
compt.human.know.raw <- read.table("../COMPARTMENTS/human_compartment_knowledge_full.tsv", 
		header = FALSE,
		sep = "\t",
		quote = "",
		comment.char = "",
		stringsAsFactors = FALSE)
colnames(compt.human.know.raw) <- compt.template.colnames
# -db_2- mouse knowledge location database
compt.mouse.know.raw <- read.table("../COMPARTMENTS/mouse_compartment_knowledge_full.tsv", 
		header = FALSE,
		sep = "\t",
		quote = "",
		comment.char = "",
		stringsAsFactors = FALSE)
colnames(compt.mouse.know.raw) <- compt.template.colnames


{   # only use once
	DoPartUnique <- function(xxpairs, cols.select=c(1:2)) 
	{
		if (sum(cols.select %in% c(1:ncol(xxpairs))) != length(cols.select)) {
			stop("Columns selected are undefined! Please check again!")
		}
		rownames(xxpairs) <- NULL
		tmp.uni <- xxpairs[, cols.select]
		tmp.uni <- unique(tmp.uni)
		xxpairs <- xxpairs[as.integer(rownames(tmp.uni)),]
		xxpairs
	}
	# check if the Ensembl Protein.stable.ID and Gene.name is same as the existing database
	compt.check.db <- compt.human.know.raw
	ensembl.check.db <- ensembl.9606.db
	# get unique list of genes in compt.*
	compt.check.sdf <- unique(compt.check.db[, c("Protein.stable.ID", "Gene.name")])
	compt.check.sdf[, c("CHECK.SRC")] <- "compt"
	# get the same 3 cols in ensembl
	ensembl.check.db <- ensembl.check.db[, c("Protein.stable.ID", "Gene.name")]
	ensembl.check.db[, c("CHECK.SRC")] <- "ensembl"
	# use DoPartUnique
	merge.check.db <- rbind(ensembl.check.db, compt.check.sdf)
	merge.check.db <- DoPartUnique(merge.check.db, c(1:2))
	# NOT the same
	# use our ensembl to 
}


# assign 
#
ensembl.tax.db <- ensembl.10090.db
entrez.tax.db <- entrez.10090.db
# process preparation
compt.proc.db <- compt.mouse.know.raw





# go.ids final mapping list
# Nucleus Cytoplasm Cytosol Cytoskeleton
# Peroxisome Lysosome Endoplasmic-reticulum Golgi-apparatus
# Plasma-membrane Endosome Extracellular-region Mitochondrion
go.targets.match <- c("GO:0005634", "GO:0005737", "GO:0005829", "GO:0005856", 
			"GO:0005777", "GO:0005764", "GO:0005783", "GO:0005794",
			"GO:0005886", "GO:0005768", "GO:0005576", "GO:0005739")
#go.ref.relation <- as.list(GO.db::GOCCANCESTOR)  # CC cellular_component, ancestor: all parent nodes

# 1st: the primary process
compt.proc.pid.fac <- levels(factor(compt.proc.db[, c("Protein.stable.ID")]))
# compt.go.ext.not.exist <- character()  [FURTURE USE]
prog.bar.slim.compt <- progress::progress_bar$new(total = length(compt.proc.pid.fac))
prog.bar.slim.compt$tick(0)
compt.go.slim.prim <- lapply(compt.proc.pid.fac, 
	compt.proc.db = compt.proc.db,
	go.cc.ancestor.ref = as.list(GO.db::GOCCANCESTOR),
	go.targets.match = go.targets.match,
	FUN = function(x, compt.proc.db, go.cc.ancestor.ref, go.targets.match) {
		this.gene.all.records <- compt.proc.db[which(compt.proc.db[, c("Protein.stable.ID")] == x), ]
		# split records by Source & Evidence
		this.rec.factor <- unique(this.gene.all.records[, c("Source", "Evidence", "score")])
		# slim the records as one evidence is put in only one row
		this.slim.records <- apply(this.rec.factor, MARGIN = 1, 
			records.all = this.gene.all.records, 
			go.cc.ancestor.ref = go.cc.ancestor.ref,
			go.targets.match = go.targets.match,
			FUN = function(x, records.all, go.cc.ancestor.ref, go.targets.match) {
				this.source <- as.character(x["Source"])
				this.evidence <- as.character(x["Evidence"])
				this.score <- as.numeric(x["score"])
				inds.src <- which(records.all[, "Source"] == this.source)
				inds.evd <- which(records.all[, "Evidence"] == this.evidence)
				this.sv.records <- records.all[intersect(inds.src, inds.evd), ]  # select one evidence each time
				this.sv.go.ids <- this.sv.records[, c("GO.ID")]
				# find 2 things, one in go.targets.match, one as the most child-like one
				this.sv.go.ancs <- lapply(this.sv.go.ids, go.cc.ancestor.ref = go.cc.ancestor.ref,
					FUN = function(x, go.cc.ancestor.ref) {
						ind.list <- which(names(go.cc.ancestor.ref) == x)
						tmp.res <- NULL
						if (length(ind.list) > 0) {
							tmp.res <- go.cc.ancestor.ref[[ind.list]]
						} else {
							warning(paste0("This GO.ID: ", x," is not exist in cur-db!"))
							# compt.go.ext.not.exist <<- c(compt.go.ext.not.exist, x)  [FURTURE USE]
						}
						tmp.res
					}
				)
				this.sv.go.ancs.uniq <- unique(as.character(unlist(this.sv.go.ancs)))
				# result - most-children ones
				this.last.children.goid <- setdiff(this.sv.go.ids, this.sv.go.ancs.uniq)  # get the most-children IDs
				if (length(this.last.children.goid) == 0) {  # it is implied that the length >= 1, so give stop exception here.
					stop(paste0("Unexpected child not exists for", records.all[1, "Gene.name"], "!"))
				}
				# result - mapping to predefined controlled vocabulary
				this.map.goid.pl <- character()
				this.last.c.goid.pl <- character()
				inds.target.m <- match(this.last.children.goid, this.sv.go.ids)
				for (i in 1:length(inds.target.m)) {
					tmp.go.ac.list <- c(this.last.children.goid[i], this.sv.go.ancs[[inds.target.m[i]]])  # fetch its ancestors and pack itself in
					tmp.tgm <- which(tmp.go.ac.list %in% go.targets.match)
					if (length(tmp.tgm) == 0) {  # not matched in predefined ctrl vocabulary, use NA
						this.map.goid.pl <- c(this.map.goid.pl, NA)
						this.last.c.goid.pl <- c(this.last.c.goid.pl, this.last.children.goid[i])
					} else {
						if (length(tmp.tgm) == 1) {
							this.map.goid.pl <- c(this.map.goid.pl, tmp.go.ac.list[tmp.tgm])
							this.last.c.goid.pl <- c(this.last.c.goid.pl, this.last.children.goid[i])
						} else {  # sometimes it maps to >1 targets
							tmp.matches <- tmp.go.ac.list[tmp.tgm]
							tmp.m.ancestor <- lapply(tmp.matches, go.cc.ancestor.ref = go.cc.ancestor.ref,
								FUN = function(x, go.cc.ancestor.ref) {
									tmp.res <- go.cc.ancestor.ref[[which(names(go.cc.ancestor.ref) == x)]]  # here no need to check as its source is selected ahead
									# add special rules
									# 1 GO:0005856 Cytoskeleton, whose ancestor doesn't include GO:0005737 Cytoplasm
									if (x == "GO:0005856") {
										tmp.res <- c(tmp.res, "GO:0005737")
									}
									tmp.res
								}
							)
							tmp.m.uniq <- unique(as.character(unlist(tmp.m.ancestor)))
							tmp.last.child <- setdiff(tmp.matches, tmp.m.uniq)
							# check
							if (sum(tmp.last.child %in% go.targets.match) != length(tmp.last.child)) {
								warning(paste0("Some of the Last Children are out of list: ", paste0(setdiff(tmp.matches, go.targets.match), collapse = ", "), "."))
							}
							if (length(tmp.last.child) > 1) {  # sometimes it exists. Map to multi-targets, see below
								#stop(paste0("Unexpected multi-children exists with ", this.sv.go.ids[inds.target.m[i]], ", and gene: ",
								#	records.all[1, "Gene.name"], ", children are: ", paste0(tmp.last.child, collapse = ", ")))
							}
							if (length(tmp.last.child) == 0) {
								warning(paste0("Unexpected zero last child, Gene: ", records.all[1, "Gene.name"]))
							}
							this.map.goid.pl <- c(this.map.goid.pl, tmp.last.child)  # use the one or multi children
							this.last.c.goid.pl <- c(this.last.c.goid.pl, rep(this.last.children.goid[i], times = length(tmp.last.child)))  # times it to match length with the above
						}
					}
				}
				# get corresponding terms
				this.map.go.term.pl <- sapply(this.map.goid.pl, function(x) {
					res <- NA
					if (!is.na(x) && !is.null(GOTERM[[x]])) {
						res <- Term(GOTERM[[x]])
					}
					res
				})
				this.map.go.term.pl <- as.character(this.map.go.term.pl)
				# get corresponding terms
				this.last.c.go.term.pl <- sapply(this.last.c.goid.pl, function(x) {
					res <- NA
					if (!is.na(x) && !is.null(GOTERM[[x]])) {
						res <- Term(GOTERM[[x]])
					}
					res
				})
				this.last.c.go.term.pl <- as.character(this.last.c.go.term.pl)
				# return
				data.frame(Protein.stable.ID = records.all[1, "Protein.stable.ID"],
					Gene.name = records.all[1, "Gene.name"],
					GO.ID.target = this.map.goid.pl,
					GO.Term.target = this.map.go.term.pl,
					GO.ID.last.child = this.last.c.goid.pl,
					GO.Term.last.child = this.last.c.go.term.pl,
					Source = this.source,
					Evidence = this.evidence,
					score = this.score,
					stringsAsFactors = FALSE
				)
			}
		)
		this.res.records <- dplyr::bind_rows(this.slim.records)
		prog.bar.slim.compt$tick()  # count once
		# return
		this.res.records
	}
)
# merge rows together
compt.go.slim.res <- dplyr::bind_rows(compt.go.slim.prim)
# get <NA>s in GOterm
inds.NAs <- which(is.na(compt.go.slim.res$GO.Term.target))
compt.go.slim.res[inds.NAs, "GO.Term.target"] <- "other"  # saveRDS
#compt.go.slim.res[, "Gene.name"] <- NULL  # remove for preparation of merging
# so, ensembl.gene name will be used, but those cannot be settled will use the name from the compt.*
compt.go.slim.f1 <- dplyr::left_join(compt.go.slim.res, ensembl.tax.db, by = c("Protein.stable.ID" = "Protein.stable.ID"))
compt.go.slim.f1[, "Gene.name"] <- compt.go.slim.f1[, "Gene.name.y"]
ind.tmp.NAs <- which(is.na(compt.go.slim.f1[, "Gene.name"]))
compt.go.slim.f1[ind.tmp.NAs, "Gene.name"] <- compt.go.slim.f1[ind.tmp.NAs, "Gene.name.x"]
# merge
compt.go.slim.f2 <- dplyr::left_join(compt.go.slim.f1, entrez.tax.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
	by = c("Gene.name" = "Symbol_from_nomenclature_authority"))
# final result
compt.go.slim.final <- compt.go.slim.f2[, c("Protein.stable.ID", "GeneID", "Gene.name",
	"GO.ID.target", "GO.Term.target", "GO.ID.last.child", "GO.Term.last.child",
	"Source", "Evidence", "score")]
# added at 2020.12.20
# by exploring the result data.
# GET: Gene.name = ENSP00000485133, or someelse like it; it is illegal and will be removed
# GET: GeneID has NAs, check the corresponding Gene.name; Decide to abandon all those with NA GeneID.
compt.go.slim.final <- compt.go.slim.final[which(!is.na(compt.go.slim.final$GeneID)), ]
# finally, get 79941 rows -> 79628 rows, human
# finally, get 66713 rows -> 66569 rows, mouse

# Other changes, added at 2020.03.16
# .simpleCaps the GO.Term.target
Tc.Cap.simple <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
}

Tc.Cap.simple.vec <- function(to.cap.vec) {
	unlist(lapply(to.cap.vec, FUN = Tc.Cap.simple))
}

compt.go.slim.final$GO.Term.target <- Tc.Cap.simple.vec(compt.go.slim.final$GO.Term.target)

# get the final result
# compt.go.slim.final
