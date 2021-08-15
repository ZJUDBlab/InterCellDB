
# This file is to compare InterCellDB with other 5 databases
library(dplyr)
library(ggplot2)
library(cowplot)
# option library
library(clusterProfiler)  # for GO enrichment analysis


# global variables
# DB are arranged by presenting time
ktDB.namelist <- c("iTALK", "CellPhoneDB", "SingleCellSignalR", "NicheNet", "CellChatDB")
ktDB.favor.color <- c("#EF98CE", "#CEC3F2", "#D9C5A4", "#96CAEA", "#8CACAF")
# InterCellDB are ordered BY # experimentally validated # from curated database # by prediction # All included
ktIT.namelist <- c("InterCellDB.exp", "InterCellDB.know", "InterCellDB.pred", "InterCellDB")
ktIT.favor.color <- c("#CDD3D4", "#DEE3E4", "#EFF3F4", "#BCC3C4")
#
ktDB.loc.all <- c("Extracellular Region", "Plasma Membrane", 
	"Cytoplasm", "Cytosol", "Cytoskeleton", 
	"Endoplasmic Reticulum", "Golgi Apparatus", 
	"Endosome", "Lysosome", "Peroxisome",
	"Mitochondrion", "Nucleus", "Other")
ktDB.loc.slim <- c("Extracellular Region", "Plasma Membrane", 
	"Cytoplasm", "Nucleus", "Other")
# actions
ktDB.act.mode <- c("activation", "binding", "catalysis", "expression", "inhibition", "ptmod", "reaction")
ktDB.act.effect <- c("positive", "negative", "unspecified", "undirected")

# picked GO terms
ktGO.picked.terms <- c("angiogenesis", "EMT", "complement activation", "cytokine activity")
	## GO:0001525 [angiogenesis] <BP>
	go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0001837 [epithelial to mesenchymal transition] <BP>
	go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0006955 [immune response]
	#go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.human.ref.db, go.human.ref.db, 
	#	go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0006956 [complement activation] <BP>
	go.comp.act.genes <- Tool.FindGenesFromGO("GO:0006956", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0005125 [cytokine activity] <MF>
	go.cytokine.act.genes <- Tool.FindGenesFromGO("GO:0005125", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
ktGO.picked.list <- list(go.vessel.genes, go.emt.genes, go.comp.act.genes, go.cytokine.act.genes)
names(ktGO.picked.list) <- ktGO.picked.terms

# backup GO terms
c("cell chemotaxis")  # <BP>
# <MF>
c("cytokine activity", "GO:0005125")
c("receptor ligand activity", "GO:0048018")
c("growth factor activity", "GO:0008083")
# <CC>
c("plasma membrane protein complex", "GO:0098797")
c("extracellular matrix", "GO:0031012")
c("receptor complex", "GO:0043235")


# - picked GO terms end - 

# tool functions for main ctrl or accessory files
#
ktgp.cut <- ">"  # split for 2 genes in gene pairs
#
Ctrl.check.inter.overlap <- function(
	cmp.list  # names must be included
) {
	this.res.list <- list()
	for (i in seq_along(cmp.list)) {  # goes from y, pick one DB-Y each time
		this.dummy.name <- names(cmp.list)[i]
		for (j in seq_along(cmp.list)) {
			tmp.cmp.name <- names(cmp.list)[j]
			tmp.inter.cnt <- length(intersect(cmp.list[[i]], cmp.list[[j]]))
			this.res.list <- c(this.res.list, list(data.frame(DBx = tmp.cmp.name, DBy = this.dummy.name, 
				inter.cnt = tmp.inter.cnt,
				inter.percent = tmp.inter.cnt / length(cmp.list[[i]]),
				stringsAsFactors = FALSE)))
		}
	}
	bind_rows(this.res.list)
} 

### Brief explanation
# [part xxx] is incremental, and never removed. Any new idea will get one new part ID.
# [deter <a-z>] is to denote finally picked ones.



##### For Describing GENE - Human
#
# 
tmp.genes.cnt.all <- list(cellphone.genes.res, nichenet.genes.res, italk.genes.res, scsr.genes.res, cellchat.human.genes.res, 
	#intercelldb.exp.human.genes.res, intercelldb.know.human.genes.res, intercelldb.pred.human.genes.res, 
	intercelldb.human.genes.res)
names(tmp.genes.cnt.all) <- c(ktDB.namelist, "InterCellDB")


# [part 001]: involved genes, count comparison
# [TODO] if need other database as reference, and compare to the rest database
	tmp.genes.count <- sapply(tmp.genes.cnt.all, function(x) {length(x)})
	tmp.genes.count.cmp.df <- data.frame(DBname = names(tmp.genes.cnt.all), 
		genes.count = tmp.genes.count, stringsAsFactors = FALSE)
	tmp.genes.count.colour <- c(ktDB.favor.color, "#BCC3C4")
	names(tmp.genes.count.colour) <- c(ktDB.namelist, "InterCellDB")
	tmp.genes.count.cmp.df$DBname <- factor(tmp.genes.count.cmp.df$DBname, 
		levels = c(ktDB.namelist, "InterCellDB"))
	pic.genes.count.cmp <- ggplot(data = tmp.genes.count.cmp.df, 
		aes(x = DBname, y = genes.count, fill = DBname)) + 
		geom_col() + 
		scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
		scale_fill_manual(values = tmp.genes.count.colour) + 
		xlab("DBname") + ylab("involved genes count") +
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
			legend.position = "bottom")


# [part 002]: gene overlaps, by refering InterCellDB (table decribing why genes not in our database should be given[TODO])
# [NOT DONE yet]

# [deprecated] [part 003]: the unique genes in InterCellDB, what they serve for
if (FALSE) {  # [deprecated] using GO enrichment analysis to illustrate
	other5.human.genes <- unique(Reduce(union, tmp.genes.cnt.all[1:5]))
	other5.human.gdf <- left_join(data.frame(gene.name = other5.human.genes, stringsAsFactors = FALSE), 
		genes.human.ref.db$gene.ncbi.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
		by = c("gene.name" = "Symbol_from_nomenclature_authority"))
	other5.human.enrich.go <- enrichGO(gene = other5.human.gdf$GeneID, OrgDb = org.Hs.eg.db, ont = "ALL", 
		pvalueCutoff = .01, qvalueCutoff = .05, pool = TRUE)
	dotplot(other5.human.enrich.go, showCategory = 20)


	uq.intercelldb.human.genes <- setdiff(intercelldb.human.genes.res, other5.human.genes)
	uq.intercelldb.human.gdf <- left_join(data.frame(gene.name = uq.intercelldb.human.genes, stringsAsFactors = FALSE), 
		genes.human.ref.db$gene.ncbi.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
		by = c("gene.name" = "Symbol_from_nomenclature_authority"))
	uq.intercelldb.human.enrich.go <- enrichGO(gene = uq.intercelldb.human.gdf$GeneID, OrgDb = org.Hs.eg.db, ont = "ALL", 
		pvalueCutoff = .01, qvalueCutoff = .05, pool = TRUE)
	dotplot(uq.intercelldb.human.enrich.go, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale = "free")
}


#
# [part 004]: check attributes of genes in other database
#
# [part 004]-sub1: subcellular location (scores selection)
	# tool function
	Ctrl.replace.locs.cytoplasm <- function(need.replace) {
		res.to.replace <- "Cytoplasm"
		if (need.replace == TRUE) {
			res.to.replace <- c(res.to.replace, c("Cytoplasm", "Cytoskeleton", "Cytosol", "Endoplasmic Reticulum", 
				"Endosome", "Endosome", "Golgi Apparatus", "Lysosome", "Mitochondrion", "Peroxisome"))
		}
		return(res.to.replace)	
	}
	#
	check.loc.each.gene.most.conf <- tapply(anno.location.human.ref.db$score, anno.location.human.ref.db$Gene.name, max)
	#   2    3    4    5   run at 2021.02.23
	# 106 3086 4482 9972
	# so we use the most confident locations for every gene by score selection
	check.loc.egm.df <- data.frame(gene = names(check.loc.each.gene.most.conf), loc.max.score = check.loc.each.gene.most.conf, stringsAsFactors = FALSE)
	check.loc.egm.df <- left_join(check.loc.egm.df, anno.location.human.ref.db[, c("GeneID", "Gene.name", "GO.Term.target", "score")],
		by = c("gene" = "Gene.name", "loc.max.score" = "score"))
	check.loc.egm.df <- unique(check.loc.egm.df)  # 28197 rows
	# 
	tmp.use.loc.category <- ktDB.loc.all  # change to use for different purpose
	tmp.use.loc.replace <- TRUE
	tmp.use.DB.namelist <- c(ktDB.namelist, "InterCellDB")
	#
	tmp.check.loc.alldb <- lapply(tmp.genes.cnt.all[tmp.use.DB.namelist], 
		ref.loc.kinds = tmp.use.loc.category, gene.loc.ref.db = check.loc.egm.df, 
		need.replace = tmp.use.loc.replace, 
		function(x, ref.loc.kinds, gene.loc.ref.db, need.replace) {
			as.integer(unlist(lapply(ref.loc.kinds, use.genes = x, gene.loc.ref.db = gene.loc.ref.db, 
				function(y, use.genes, gene.loc.ref.db) {
					if (y == "Cytoplasm") {
						y <- Ctrl.replace.locs.cytoplasm(need.replace = need.replace)
					}
					this.ref.genes <- gene.loc.ref.db[which(gene.loc.ref.db$GO.Term.target %in% y), "gene"]
					length(intersect(use.genes, this.ref.genes))
				})))
		})
	tmp.loc.anal.df <- data.frame(DBname = rep(tmp.use.DB.namelist, each = length(tmp.use.loc.category)),
		Location = rep(tmp.use.loc.category, times = length(tmp.use.DB.namelist)),
		count = Reduce(c, tmp.check.loc.alldb),
		stringsAsFactors = FALSE)
	tmp.loc.anal.df$count <- sapply(tmp.loc.anal.df$count, function(x) {if (x > 1000) x <- 1000; x}, USE.NAMES = FALSE)
	# remove InterCellDB
	#tmp.loc.anal.df <- tmp.loc.anal.df[which(tmp.loc.anal.df$DBname != "InterCellDB"), ]
	# add factor
	tmp.loc.anal.df$Location <- factor(tmp.loc.anal.df$Location, levels = tmp.use.loc.category)
	tmp.loc.anal.df$DBname <- factor(tmp.loc.anal.df$DBname, levels = tmp.use.DB.namelist)
	tmp.loc.anal.plot <- ggplot(data = tmp.loc.anal.df) + 
			geom_raster(aes(x = DBname, y = Location, fill = count)) + 
			scale_fill_gradient(name = "Count", low = "yellow", high = "red") + 
			theme_cowplot(16) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# [part 004]-sub2: molecular function from Uniprot
	tmp.use.mf.unip.sel <- c("Cytokine", "Growth Factor", "Hormone", "Ion Channel", "Transducer")
	# use Uniprot remap to mapping those keys to get its subset of Keywords
	tmp.use.mf.unip.genelist <- lapply(tmp.use.mf.unip.sel, 
		ref.unip.remap = Uniprot.key.map.list,
		ref.unip.anno.db = anno.type.human.ref.db, 
		function(x, ref.unip.remap, ref.unip.anno.db) {
			this.unip.keys <- ref.unip.remap[which(ref.unip.remap[, 2] == x), 1]
			unique(ref.unip.anno.db[which(ref.unip.anno.db$Keyword.Name %in% this.unip.keys), "Gene.name"])
		})
	names(tmp.use.mf.unip.genelist) <- tmp.use.mf.unip.sel
	# run Uniprot Keywords related genes and DB genes intersection
	tmp.use.DB.namelist.for.unip <- c(ktDB.namelist, "InterCellDB")
	tmp.DBgenes.unip.countlist <- lapply(seq_along(tmp.genes.cnt.all[tmp.use.DB.namelist.for.unip]),
		inside.genes.cnt.all = tmp.genes.cnt.all[tmp.use.DB.namelist.for.unip], unip.list = tmp.use.mf.unip.genelist,
		function(x, inside.genes.cnt.all, unip.list) {
			this.db.genes <- inside.genes.cnt.all[[x]]
			this.db.unip.inter <- vapply(unip.list, it.genes = this.db.genes,
				FUN.VALUE = integer(1),
				USE.NAMES = FALSE,
				FUN = function(y, it.genes) {
					length(intersect(y, it.genes))
					})
			data.frame(DBname = rep(names(inside.genes.cnt.all)[x], times = length(unip.list)),
				Unip.Keywords = names(unip.list),
				unip.cnt = this.db.unip.inter,
				stringsAsFactors = FALSE)
		})
	tmp.DBgenes.unip.countdf <- bind_rows(tmp.DBgenes.unip.countlist)
	tmp.DBgenes.unip.countdf$DBname <- factor(tmp.DBgenes.unip.countdf$DBname, levels = tmp.use.DB.namelist.for.unip)
	tmp.for.unip.colour <- c(ktDB.favor.color, ktIT.favor.color)
	names(tmp.for.unip.colour) <- c(ktDB.namelist, "InterCellDB")
	# plot
	tmp.for.unip.plot <- ggplot(data = tmp.DBgenes.unip.countdf) + 
		geom_col(aes(x = DBname, y = unip.cnt, group = Unip.Keywords, fill = DBname)) + 
		facet_grid(.~Unip.Keywords) + 
		scale_fill_manual(values = tmp.for.unip.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			legend.position = "bottom")


# [part 004]-sub3: picked GO terms
	tmp.chk.genes.picked.GOlist <- list(
		# <MF>
		c("cytokine activity", "GO:0005125"),
		c("receptor ligand activity", "GO:0048018"),
		#c("signaling receptor activity", "GO:0038023"),
		#c("growth factor activity", "GO:0008083"),
		# <CC>
		c("receptor complex", "GO:0043235"),
		c("plasma membrane protein complex", "GO:0098797"),
		c("extracellular matrix", "GO:0031012"),
		# <BP>
		c("angiogenesis", "GO:0001525"),
		c("EMT", "GO:0001837"),  # epithelial to mesenchymal transition
		c("complement activation", "GO:0006956"),
		c("cell chemotaxis", "GO:0060326")
	)
	# subset the GOlist
	tmp.chk.genes.picked.GOlist <- tmp.chk.genes.picked.GOlist[6:9]
	#
	# get genes
	tmp.chk.genes.picked.goid <- vapply(tmp.chk.genes.picked.GOlist, function(x) {x[2]}, FUN.VALUE = character(1))
	tmp.picked.go.rel.genelist <- lapply(tmp.chk.genes.picked.goid, 
		function(x) {
			Tool.FindGenesFromGO(x, genes.human.ref.db, go.human.ref.db, 
				go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
		})
	names(tmp.picked.go.rel.genelist) <- vapply(tmp.chk.genes.picked.GOlist, function(x) {x[1]}, FUN.VALUE = character(1))

	# run GO related genes and DB genes intersection
	tmp.use.DB.namelist.picked.GO <- c(ktDB.namelist, "InterCellDB")
	tmp.DBgenes.GO.countlist <- lapply(seq_along(tmp.genes.cnt.all[tmp.use.DB.namelist.picked.GO]), 
		inside.genes.cnt.all = tmp.genes.cnt.all[tmp.use.DB.namelist.picked.GO], GO.list = tmp.picked.go.rel.genelist, 
		function(x, inside.genes.cnt.all, GO.list) {
			this.db.genes <- inside.genes.cnt.all[[x]]
			this.db.go.inter <- vapply(GO.list, it.genes = this.db.genes, 
				FUN.VALUE = integer(1), 
				USE.NAMES = FALSE, 
				FUN = function(y, it.genes) {
					length(intersect(y, it.genes))
				})
			data.frame(DBname = rep(names(inside.genes.cnt.all)[x], times = length(GO.list)),
				GOterm = names(GO.list),
				go.cnt = this.db.go.inter,
				stringsAsFactors = FALSE)
		})
	tmp.DBgenes.GO.countdf <- bind_rows(tmp.DBgenes.GO.countlist)
	tmp.DBgenes.GO.countdf$DBname <- factor(tmp.DBgenes.GO.countdf$DBname, levels = tmp.use.DB.namelist.picked.GO)
	tmp.DBgenes.GO.countdf$GOterm <- factor(tmp.DBgenes.GO.countdf$GOterm, levels = names(tmp.picked.go.rel.genelist))
	tmp.for.GO.colour <- c(ktDB.favor.color, ktIT.favor.color)
	names(tmp.for.GO.colour) <- c(ktDB.namelist, "InterCellDB")
	# plot
	tmp.go.explore.plot <- ggplot(data = tmp.DBgenes.GO.countdf) + 
		geom_col(aes(x = DBname, y = go.cnt, group = GOterm, fill = DBname)) + 
		facet_grid(.~GOterm) + 
		scale_fill_manual(values = tmp.for.GO.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			strip.text.x = element_text(size = 6), 
			legend.position = "bottom")
	# width = 8.22, height = 5.65

#
# [part 005]: check ligand and receptor overlaps for databases
# [part 005]-sub1: chk receptor
tmp.receptor.cnt.all <- list(tmp.receptor.italk.res, tmp.receptor.cpdb.res, tmp.receptor.scsr.res, 
	tmp.receptor.nichenet.res, tmp.receptor.cellchat.res, tmp.receptor.intercelldb.res)
names(tmp.receptor.cnt.all) <- c(ktDB.namelist, "InterCellDB")
#
tmp.inter.receptor.df <- Ctrl.check.inter.overlap(tmp.receptor.cnt.all)
# prep plot
tmp.inter.receptor.df$DBx <- factor(tmp.inter.receptor.df$DBx, levels = names(tmp.receptor.cnt.all))
tmp.inter.receptor.df$DBy <- factor(tmp.inter.receptor.df$DBy, levels = names(tmp.receptor.cnt.all))
tmp.inter.receptor.df$inter.percent <- round(tmp.inter.receptor.df$inter.percent, digits = 2)
tmp.inter.receptor.df$inter.percent.label <- format(tmp.inter.receptor.df$inter.percent, scientific = FALSE)
# plot
tmp.rec.inter.plot <- ggplot(data = tmp.inter.receptor.df) + 
	geom_raster(aes(x = DBx, y = DBy, fill = inter.percent)) + 
	geom_text(aes(x = DBx, y = DBy, label = inter.percent.label)) + 
	scale_fill_gradient(low = "yellow", high = "red") +
	theme_cowplot(16) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
# (width = 6.5, height=5.4)


# [part 005]-sub2: chk ligand
tmp.ligand.cnt.all <- list(tmp.ligand.italk.res, tmp.ligand.cpdb.res, tmp.ligand.scsr.res, 
	tmp.ligand.nichenet.res, tmp.ligand.cellchat.res, tmp.ligand.intercelldb.res)
names(tmp.ligand.cnt.all) <- c(ktDB.namelist, "InterCellDB")
#
tmp.inter.ligand.df <- Ctrl.check.inter.overlap(tmp.ligand.cnt.all)
# prep plot
tmp.inter.ligand.df$DBx <- factor(tmp.inter.ligand.df$DBx, levels = names(tmp.ligand.cnt.all))
tmp.inter.ligand.df$DBy <- factor(tmp.inter.ligand.df$DBy, levels = names(tmp.ligand.cnt.all))
tmp.inter.ligand.df$inter.percent <- round(tmp.inter.ligand.df$inter.percent, digits = 2)
tmp.inter.ligand.df$inter.percent.label <- format(tmp.inter.ligand.df$inter.percent, scientific = FALSE)
# plot
tmp.rec.inter.plot <- ggplot(data = tmp.inter.ligand.df) + 
	geom_raster(aes(x = DBx, y = DBy, fill = inter.percent)) + 
	geom_text(aes(x = DBx, y = DBy, label = inter.percent.label)) + 
	scale_fill_gradient(low = "yellow", high = "red") +
	theme_cowplot(16) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
# (width = 6.5, height=5.4)



# ------ preserve some part IDs for extension in GENE ------


##### For Describing GENE PAIRS - Human
#
#
cmp.pairs.db.list <- list(remap.italk.pairs.human, remap.cellphone.pairs.human, remap.scsr.pairs.human, 
	remap.nichenet.pairs.human, remap.cellchat.pairs.human,
	ref.intercelldb.pairs.human.experiments, ref.intercelldb.pairs.human.knowledge, 
	ref.intercelldb.pairs.human.prediction, ref.intercelldb.pairs.human)
	#ref.intercelldb.pairs.human.LR)
names(cmp.pairs.db.list) <- c(ktDB.namelist, ktIT.namelist)


# [part 020]: gene pairs total count comparison
	gpairs.count.cmp.df <- bind_rows(lapply(seq_along(cmp.pairs.db.list), use.tg = cmp.pairs.db.list, 
		function(x, use.tg) {
			data.frame(DBname = names(use.tg)[x], gpairs.cnt = length(use.tg[[x]]))
		}))
	gpairs.count.cmp.df$DBname <- factor(gpairs.count.cmp.df$DBname, levels = names(cmp.pairs.db.list))
	#
	gpairs.count.cmp.df$gpairs.sqrt.cnt <- sqrt(gpairs.count.cmp.df$gpairs.cnt)
	# plot
	tmp.gpairs.count.colour <- c(ktDB.favor.color, ktIT.favor.color)
	names(tmp.gpairs.count.colour) <- c(ktDB.namelist, ktIT.namelist)
	tmp.gpairs.count.breaks <- c(0, 100, 400, 900, 1600, 2500)
	tmp.gpairs.count.plot <- ggplot(data = gpairs.count.cmp.df) + 
		geom_col(aes(x = DBname, y = gpairs.sqrt.cnt, fill = DBname)) + 
		scale_y_continuous(limits = c(0, 2500),
			breaks = tmp.gpairs.count.breaks, 
			labels = tmp.gpairs.count.breaks * tmp.gpairs.count.breaks,
			expand = expansion(mult = c(0, .02))) + 
		scale_fill_manual(values = tmp.gpairs.count.colour) + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank(), 
				legend.position = "right")
	# (width = 10, height = 5.5) test1


#  [deprecated] [part 021]: gene pairs overlap [InterCellDB being used as benchmark]
if (FALSE) {
	tmp.chk.overlap.sub0.cmp.db <- cmp.pairs.db.list[ktDB.namelist]
	tmp.chk.overlap.sub0.df <- bind_rows(lapply(seq_along(tmp.chk.overlap.sub0.cmp.db),
		to.cmp.list = tmp.chk.overlap.sub0.cmp.db, gpairs.ref = cmp.pairs.db.list, 
		function(x, to.cmp.list, gpairs.ref) {
			Ctrl.percent.plot.add.data(to.cmp.list[[x]], names(to.cmp.list)[x],
				col.desc.vals = c("match.in.exp", "match.in.know", "match.in.pred", "unmatch.all"),
				ref.vector = gpairs.ref[setdiff(ktIT.namelist, "InterCellDB")])
		}))
	pic.chk.overlap.sub0 <- ggplot(data = tmp.chk.overlap.sub0.df, 
			aes(x = target.names, y = result.cnt, fill = cmp.types)) + 
			geom_col(position = position_stack(reverse = TRUE)) +  
			scale_x_discrete(limits = ktDB.namelist, breaks = ktDB.namelist) + 
			scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
			scale_fill_brewer(palette = 1, direction = -1) + 
			xlab("Source") + ylab("Matches") +
			theme_cowplot(12)
}


# [part 022]: create something like relation matrix plot
# plot design
# ------------------------------
#        Formula:    intersect(DB-X, DB-Y) / DB-Y, How TEST matched in REF (%)
# <REF> R1, R2 
# <TEST> T1, T2
# ------------------------------
#  R1
#  R2
#  T1
#  T2   
#
#   x   R1  R2  T1  T2
# SO
#    (xR, yT) means how T get precision, means True Positive / (TP + FP)
#    (xT, yR) means how T get recall, means TP / (TP + FN)
# ====== In InterCellDB settings ======
#  treat other 5 database each as reference
#
	tmp.base.cmp.pairs.db.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.LR))
	names(tmp.base.cmp.pairs.db.list) <- c(ktDB.namelist, ktIT.namelist, "InterCellDB.LR")
	#
	tmp.gpairs.corel.list <- list()
	for (i in seq_along(tmp.base.cmp.pairs.db.list)) {  # goes from y, pick one DB-Y each time
		this.dummy.name <- names(tmp.base.cmp.pairs.db.list)[i]
		for (j in seq_along(tmp.base.cmp.pairs.db.list)) {
			tmp.cmp.name <- names(tmp.base.cmp.pairs.db.list)[j]
			tmp.inter.cnt <- length(intersect(tmp.base.cmp.pairs.db.list[[i]], tmp.base.cmp.pairs.db.list[[j]]))
			tmp.gpairs.corel.list <- c(tmp.gpairs.corel.list, list(data.frame(DBx = tmp.cmp.name, DBy = this.dummy.name, 
				inter.cnt = tmp.inter.cnt,
				inter.percent = tmp.inter.cnt / length(tmp.base.cmp.pairs.db.list[[i]]),
				stringsAsFactors = FALSE)))
		}
	}
	tmp.gpairs.corel.df <- bind_rows(tmp.gpairs.corel.list)
	tmp.gpairs.corel.df$DBx <- factor(tmp.gpairs.corel.df$DBx, levels = names(tmp.base.cmp.pairs.db.list))
	tmp.gpairs.corel.df$DBy <- factor(tmp.gpairs.corel.df$DBy, levels = names(tmp.base.cmp.pairs.db.list))
	tmp.gpairs.corel.df$inter.percent <- round(tmp.gpairs.corel.df$inter.percent, digits = 4)
	tmp.gpairs.corel.df$inter.percent.label <- format(tmp.gpairs.corel.df$inter.percent, scientific = FALSE)
	# plot
	ggplot(data = tmp.gpairs.corel.df) + 
		geom_raster(aes(x = DBx, y = DBy, fill = inter.percent)) + 
		geom_text(aes(x = DBx, y = DBy, label = inter.percent.label)) + 
		scale_fill_gradient(low = "yellow", high = "red") +
		theme_cowplot(16) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# [part 022]-sub-supple1: get score leveled pairs db
tmp.use.cmp.pairs.db.list <- list(remap.italk.pairs.human, remap.cellphone.pairs.human, remap.scsr.pairs.human, 
	remap.nichenet.pairs.human, remap.cellchat.pairs.human,
	highconf.intercelldb.pairs.human.experiments, midconf.intercelldb.pairs.human.experiments, lowconf.intercelldb.pairs.human.experiments, 
	ref.intercelldb.pairs.human.experiments,
	highconf.intercelldb.pairs.human.knowledge, midconf.intercelldb.pairs.human.knowledge, lowconf.intercelldb.pairs.human.knowledge, 
	ref.intercelldb.pairs.human.knowledge,
	highconf.intercelldb.pairs.human.prediction, midconf.intercelldb.pairs.human.prediction, lowconf.intercelldb.pairs.human.prediction,
	ref.intercelldb.pairs.human.prediction,
	ref.intercelldb.pairs.human)
names(tmp.use.cmp.pairs.db.list) <- c(ktDB.namelist, 
	paste(ktIT.namelist[1], c("high", "mid", "low"), sep = "."), ktIT.namelist[1],
	paste(ktIT.namelist[2], c("high", "mid", "low"), sep = "."), ktIT.namelist[2],
	paste(ktIT.namelist[3], c("high", "mid", "low"), sep = "."), ktIT.namelist[3],
	ktIT.namelist[4])
tmp.gpairs.supple.corel.list <- list()
for (i in seq_along(tmp.use.cmp.pairs.db.list)) {  # goes from y, pick one DB-Y each time
	this.dummy.name <- names(tmp.use.cmp.pairs.db.list)[i]
	for (j in seq_along(tmp.use.cmp.pairs.db.list)) {
		tmp.cmp.name <- names(tmp.use.cmp.pairs.db.list)[j]
		tmp.inter.cnt <- length(intersect(tmp.use.cmp.pairs.db.list[[i]], tmp.use.cmp.pairs.db.list[[j]]))
		tmp.gpairs.supple.corel.list <- c(tmp.gpairs.supple.corel.list, list(data.frame(DBx = tmp.cmp.name, DBy = this.dummy.name, 
			inter.cnt = tmp.inter.cnt,
			inter.percent = tmp.inter.cnt / length(tmp.use.cmp.pairs.db.list[[i]]),
			stringsAsFactors = FALSE)))
	}
}
tmp.gpairs.supple.corel.df <- bind_rows(tmp.gpairs.supple.corel.list)
tmp.gpairs.supple.corel.df$DBx <- factor(tmp.gpairs.supple.corel.df$DBx, levels = names(tmp.use.cmp.pairs.db.list))
tmp.gpairs.supple.corel.df$DBy <- factor(tmp.gpairs.supple.corel.df$DBy, levels = names(tmp.use.cmp.pairs.db.list))
tmp.gpairs.supple.corel.df$inter.percent <- round(tmp.gpairs.supple.corel.df$inter.percent, digits = 3)
tmp.gpairs.supple.corel.df$inter.percent.label <- format(tmp.gpairs.supple.corel.df$inter.percent, scientific = FALSE)
# plot
ggplot(data = tmp.gpairs.supple.corel.df) + 
	geom_raster(aes(x = DBx, y = DBy, fill = inter.percent)) + 
	geom_text(aes(x = DBx, y = DBy, label = inter.percent.label)) + 
	scale_fill_gradient(low = "yellow", high = "red") +
	theme_cowplot(16) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# [part 023]: annotate other DBs with action mode and action effects
# Q: one gene pair could have multiple modes as well as effects?
# A: For Mode, Put activation and inhibition in the main plot, others given in supple
#    For Effect, use all 4 categories, add gives 0-1 percentage
#
# [To Note] as action mode and action effect get many overlaps, so the matching process below is 
# firstly matched by overall merged unique pairs, and then in matched count, every category gets its
# percentage by each matched count. Do so, to make sure the unmatched percentage ~stable~.
# 
tmp.use.cmp.pairs.db.list <- c(cmp.pairs.db.list[1:5], "InterCellDB" = list(ref.intercelldb.pairs.human))

# action mode
	tmp.act.mode.inter.list <- list()
	for (i in seq_along(tmp.use.cmp.pairs.db.list)) {  # select other 5 database
		this.chk.db.name <- names(tmp.use.cmp.pairs.db.list)[i]
		for (j in seq_along(ref.intercelldb.act.mode.pairs.human)) {
			tmp.count <- length(intersect(tmp.use.cmp.pairs.db.list[[i]], ref.intercelldb.act.mode.pairs.human[[j]]))
			tmp.act.mode.inter.list <- c(tmp.act.mode.inter.list, list(data.frame(DBname = this.chk.db.name,
				act.mode = names(ref.intercelldb.act.mode.pairs.human)[j],
				inter.raw.cnt  = tmp.count,
				stringsAsFactors = FALSE)))
		}
		tmp.count <- length(setdiff(tmp.use.cmp.pairs.db.list[[i]], unlist(ref.intercelldb.act.mode.pairs.human)))
		tmp.act.mode.inter.list <- c(tmp.act.mode.inter.list, list(data.frame(DBname = this.chk.db.name,
				act.mode = "UNMATCHED",
				inter.raw.cnt  = tmp.count,
				stringsAsFactors = FALSE)))
	}
	tmp.act.mode.inter.df <- bind_rows(tmp.act.mode.inter.list)
	# get recalc percentage by prioritize UNMATCHED
	tmp.percent.result <- numeric()
	for (i in seq_along(tmp.use.cmp.pairs.db.list)) {
		this.chk.db.name <- names(tmp.use.cmp.pairs.db.list)[i]
		this.chk.cnt <- length(tmp.use.cmp.pairs.db.list[[i]])
		tmp.chk.part <- tmp.act.mode.inter.df[which(tmp.act.mode.inter.df$DBname == this.chk.db.name), ]
		tmp.cnt.unmatched <- tmp.chk.part[which(tmp.chk.part$act.mode == "UNMATCHED"), "inter.raw.cnt"]
		tmp.cnt.allother <- sum(tmp.chk.part$inter.raw.cnt) - tmp.cnt.unmatched
		tmp.chk.part$inter.percent <- tmp.chk.part$inter.raw.cnt / tmp.cnt.allother * 
			((this.chk.cnt - tmp.cnt.unmatched) / this.chk.cnt)
		tmp.chk.part[which(tmp.chk.part$act.mode == "UNMATCHED"), "inter.percent"] <- tmp.cnt.unmatched / this.chk.cnt
		tmp.percent.result <- c(tmp.percent.result, tmp.chk.part$inter.percent)
	}
	tmp.act.mode.inter.df$inter.percent <- tmp.percent.result
	tmp.replace <- tmp.act.mode.inter.df$act.mode
	tmp.replace[which(tmp.replace == "UNMATCHED")] <- "undefined"
	tmp.act.mode.inter.df$act.mode <- tmp.replace
	# plot
	tmp.act.mode.fill.color <- c("#D70051", "#00913A", "#1296D4", "#956134", "#C8DC32", "#B5B5B6", "#0A0AFF", "grey")
	names(tmp.act.mode.fill.color) <- c("activation", "inhibition", "binding", "catalysis", "reaction", "expression", "ptmod", "undefined")
	ggplot(data = tmp.act.mode.inter.df, 
				aes(x = DBname, y = inter.percent, fill = act.mode)) + 
				geom_col(position = position_stack(reverse = TRUE)) +  
				scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
				scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
				scale_fill_manual(values = tmp.act.mode.fill.color) + 
				xlab("DBs") + ylab("Act.Mode.matching.percentage") +
				theme_cowplot(12)
	# plot with undefined removed
	ggplot(data = tmp.act.mode.inter.df[which(tmp.act.mode.inter.df$act.mode != "undefined"),], 
				aes(x = DBname, y = inter.percent, fill = act.mode)) + 
				geom_col(position = position_fill(reverse = TRUE)) +  
				scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
				scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
				scale_fill_manual(values = tmp.act.mode.fill.color[1:7]) + 
				xlab("DBs") + ylab("Act.Mode.matching.percentage") +
				theme_cowplot(12)

	# only show the activation and inhibition (NOT THAT GOOD) (InterCellDB do get some bit more inhibition proportion)
	ggplot(data = tmp.act.mode.inter.df[which(tmp.act.mode.inter.df$act.mode %in% c("activation", "inhibition")), ], 
			aes(x = DBname, y = inter.percent, fill = act.mode)) + 
			geom_col(position = position_fill(reverse = TRUE)) +  
			scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
			scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
			scale_fill_manual(values = tmp.act.mode.fill.color[1:7]) + 
			xlab("DBs") + ylab("Act.Mode.matching.percentage") +
			theme_cowplot(12)

# action effect
	tmp.act.effect.inter.list <- list()
	for (i in seq_along(tmp.use.cmp.pairs.db.list)) {  # select other 5 database
		this.chk.db.name <- names(tmp.use.cmp.pairs.db.list)[i]
		for (j in seq_along(ref.intercelldb.act.effect.pairs.human)) {
			tmp.count <- length(intersect(tmp.use.cmp.pairs.db.list[[i]], ref.intercelldb.act.effect.pairs.human[[j]]))
			tmp.act.effect.inter.list <- c(tmp.act.effect.inter.list, list(data.frame(DBname = this.chk.db.name,
				act.effect = names(ref.intercelldb.act.effect.pairs.human)[j],
				inter.raw.cnt  = tmp.count,
				stringsAsFactors = FALSE)))
		}
		tmp.count <- length(setdiff(tmp.use.cmp.pairs.db.list[[i]], unlist(ref.intercelldb.act.effect.pairs.human)))
		tmp.act.effect.inter.list <- c(tmp.act.effect.inter.list, list(data.frame(DBname = this.chk.db.name,
				act.effect = "UNMATCHED",
				inter.raw.cnt  = tmp.count,
				stringsAsFactors = FALSE)))
	}
	tmp.act.effect.inter.df <- bind_rows(tmp.act.effect.inter.list)
	# get recalc percentage by prioritize UNMATCHED
	tmp.percent.result <- numeric()
	for (i in seq_along(tmp.use.cmp.pairs.db.list)) {
		this.chk.db.name <- names(tmp.use.cmp.pairs.db.list)[i]
		this.chk.cnt <- length(tmp.use.cmp.pairs.db.list[[i]])
		tmp.chk.part <- tmp.act.effect.inter.df[which(tmp.act.effect.inter.df$DBname == this.chk.db.name), ]
		tmp.cnt.unmatched <- tmp.chk.part[which(tmp.chk.part$act.effect == "UNMATCHED"), "inter.raw.cnt"]
		tmp.cnt.allother <- sum(tmp.chk.part$inter.raw.cnt) - tmp.cnt.unmatched
		tmp.chk.part$inter.percent <- tmp.chk.part$inter.raw.cnt / tmp.cnt.allother * 
			((this.chk.cnt - tmp.cnt.unmatched) / this.chk.cnt)
		tmp.chk.part[which(tmp.chk.part$act.effect == "UNMATCHED"), "inter.percent"] <- tmp.cnt.unmatched / this.chk.cnt
		tmp.percent.result <- c(tmp.percent.result, tmp.chk.part$inter.percent)
	}
	tmp.act.effect.inter.df$inter.percent <- tmp.percent.result
	tmp.replace <- tmp.act.effect.inter.df$act.effect
	tmp.replace[which(tmp.replace == "UNMATCHED")] <- "undefined"
	tmp.act.effect.inter.df$act.effect <- tmp.replace
	# factor act.effect
	tmp.act.effect.inter.df$act.effect <- factor(tmp.act.effect.inter.df$act.effect, levels = c(ktDB.act.effect, "undefined"))
	# plot
	tmp.act.effect.fill.color <- c("red", "green", "yellow", "blue", "grey")
	names(tmp.act.effect.fill.color) <- c(ktDB.act.effect, "undefined")
	ggplot(data = tmp.act.effect.inter.df, 
				aes(x = DBname, y = inter.percent, fill = act.effect)) + 
				geom_col(position = position_stack(reverse = TRUE)) +  
				scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
				scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
				scale_fill_manual(values = tmp.act.effect.fill.color) + 
				xlab("DBs") + ylab("act.effect.matching.percentage") +
				theme_cowplot(12)
	# plot with undefined removed
	ggplot(data = tmp.act.effect.inter.df[which(tmp.act.effect.inter.df$act.effect != "undefined"),], 
				aes(x = DBname, y = inter.percent, fill = act.effect)) + 
				geom_col(position = position_fill(reverse = TRUE)) +  
				scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
				scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
				scale_fill_manual(values = tmp.act.effect.fill.color[1:4]) + 
				xlab("DBs") + ylab("Act.effect.matching.percentage") +
				theme_cowplot(12)
	# plot with undefined removed
	ggplot(data = tmp.act.effect.inter.df[which(tmp.act.effect.inter.df$act.effect %in% c("positive", "negative")),], 
				aes(x = DBname, y = inter.percent, fill = act.effect)) + 
				geom_col(position = position_fill(reverse = TRUE)) +  
				scale_x_discrete(limits = names(tmp.use.cmp.pairs.db.list), breaks = names(tmp.use.cmp.pairs.db.list)) + 
				scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
				scale_fill_manual(values = tmp.act.effect.fill.color[1:4]) + 
				xlab("DBs") + ylab("Act.effect.matching.percentage") +
				theme_cowplot(12)

