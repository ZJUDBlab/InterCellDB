
# This file shares plenty of settings given in Main-v2-human-process.R

# global variables
ktDB.mouse.namelist <- c("SingleCellSignalR", "NicheNet", "CellChatDB")
ktDB.mouse.favor.color <- c("#D9C5A4", "#96CAEA", "#8CACAF")
ktIT.mouse.namelist <- c("InterCellDB.exp", "InterCellDB.know", "InterCellDB.pred", "InterCellDB")
ktIT.mouse.favor.color <- c("#BCC3C4", "#BCC3C4", "#BCC3C4", "#BCC3C4")
#
ktDB.mouse.loc.all <- c("Extracellular Region", "Plasma Membrane", 
	"Cytoplasm", "Cytosol", "Cytoskeleton", 
	"Endoplasmic Reticulum", "Golgi Apparatus", 
	"Endosome", "Lysosome", "Peroxisome",
	"Mitochondrion", "Nucleus", "Other")
ktDB.mouse.loc.slim <- c("Extracellular Region", "Plasma Membrane", 
	"Cytoplasm", "Nucleus", "Other")
# actions
ktDB.mouse.act.mode <- c("activation", "binding", "catalysis", "expression", "inhibition", "ptmod", "reaction")
ktDB.mouse.act.effect <- c("positive", "negative", "unspecified", "undirected")


##### For Describing GENE - Mouse
#
# 
mmp.genes.cnt.all <- list(scsr.mouse.genes.res, nichenet.mouse.genes.res, cellchat.mouse.genes.res, 
	#intercelldb.exp.mouse.genes.res, intercelldb.know.mouse.genes.res, intercelldb.pred.mouse.genes.res, 
	intercelldb.mouse.genes.res)
names(mmp.genes.cnt.all) <- c(ktDB.mouse.namelist, "InterCellDB")


# [part 001]: involved genes, count comparison
# [TODO] if need other database as reference, and compare to the rest database
	mmp.genes.count <- sapply(mmp.genes.cnt.all, function(x) {length(x)})
	mmp.genes.count.cmp.df <- data.frame(DBname = names(mmp.genes.cnt.all), 
		genes.count = mmp.genes.count, stringsAsFactors = FALSE)
	mmp.genes.count.colour <- c(ktDB.mosue.favor.color, "#BCC3C4")
	names(mmp.genes.count.colour) <- c(ktDB.mouse.namelist, "InterCellDB")
	mmp.genes.count.cmp.df$DBname <- factor(mmp.genes.count.cmp.df$DBname, 
		levels = c(ktDB.mouse.namelist, "InterCellDB"))
	pic.mouse.genes.count.cmp <- ggplot(data = mmp.genes.count.cmp.df, 
		aes(x = DBname, y = genes.count, fill = DBname)) + 
		geom_col() + 
		scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
		scale_fill_manual(values = mmp.genes.count.colour) + 
		xlab("DBname") + ylab("involved genes count") +
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
			legend.position = "right")
# (width = 6.7, height = 5.5)

# remove [part 002] [part 003] to align with human part

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
	check.loc.mouse.each.gene.most.conf <- tapply(anno.location.mouse.ref.db$score, anno.location.mouse.ref.db$Gene.name, max)
	#  2    3    4    5    run at 2021.03.03
	# 90 2490 5194 7792 
	# so we use the most confident locations for every gene by score selection
	check.loc.mouse.egm.df <- data.frame(gene = names(check.loc.mouse.each.gene.most.conf), loc.max.score = check.loc.mouse.each.gene.most.conf, stringsAsFactors = FALSE)
	check.loc.mouse.egm.df <- left_join(check.loc.mouse.egm.df, anno.location.mouse.ref.db[, c("GeneID", "Gene.name", "GO.Term.target", "score")],
		by = c("gene" = "Gene.name", "loc.max.score" = "score"))
	check.loc.mouse.egm.df <- unique(check.loc.mouse.egm.df)  # 26082 rows
	# 
	mmp.use.loc.category <- ktDB.mouse.loc.all  # change to use for different purpose
	mmp.use.loc.replace <- TRUE
	mmp.use.DB.namelist <- c(ktDB.mouse.namelist, "InterCellDB")
	#
	mmp.check.loc.alldb <- lapply(mmp.genes.cnt.all[mmp.use.DB.namelist], 
		ref.loc.kinds = mmp.use.loc.category, gene.loc.ref.db = check.loc.mouse.egm.df, 
		need.replace = mmp.use.loc.replace, 
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
	mmp.loc.anal.df <- data.frame(DBname = rep(mmp.use.DB.namelist, each = length(mmp.use.loc.category)),
		Location = rep(mmp.use.loc.category, times = length(mmp.use.DB.namelist)),
		count = Reduce(c, mmp.check.loc.alldb),
		stringsAsFactors = FALSE)
	mmp.loc.anal.df$count <- sapply(mmp.loc.anal.df$count, function(x) {if (x > 1000) x <- 1000; x}, USE.NAMES = FALSE)
	# remove InterCellDB
	#mmp.loc.anal.df <- mmp.loc.anal.df[which(mmp.loc.anal.df$DBname != "InterCellDB"), ]
	# add factor
	mmp.loc.anal.df$Location <- factor(mmp.loc.anal.df$Location, levels = mmp.use.loc.category)
	mmp.loc.anal.df$DBname <- factor(mmp.loc.anal.df$DBname, levels = mmp.use.DB.namelist)
	mmp.loc.anal.plot <- ggplot(data = mmp.loc.anal.df) + 
			geom_raster(aes(x = DBname, y = Location, fill = count)) + 
			scale_fill_gradient(name = "count", low = "yellow", high = "red") + 
			theme_cowplot(16) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	# (width = 5.9, height = 5.9)

# [part 004]-sub2: molecular function from Uniprot
	mmp.use.mf.unip.sel <- c("Cytokine", "Growth Factor", "Hormone", "Ion Channel", "Transducer")
	# use Uniprot remap to mapping those keys to get its subset of Keywords
	mmp.use.mf.unip.genelist <- lapply(mmp.use.mf.unip.sel, 
		ref.unip.remap = Uniprot.key.map.list,
		ref.unip.anno.db = anno.type.mouse.ref.db, 
		function(x, ref.unip.remap, ref.unip.anno.db) {
			this.unip.keys <- ref.unip.remap[which(ref.unip.remap[, 2] == x), 1]
			unique(ref.unip.anno.db[which(ref.unip.anno.db$Keyword.Name %in% this.unip.keys), "Gene.name"])
		})
	names(mmp.use.mf.unip.genelist) <- mmp.use.mf.unip.sel
	# run Uniprot Keywords related genes and DB genes intersection
	mmp.use.DB.namelist.for.unip <- c(ktDB.mouse.namelist, "InterCellDB")
	mmp.DBgenes.unip.countlist <- lapply(seq_along(mmp.genes.cnt.all[mmp.use.DB.namelist.for.unip]),
		inside.genes.cnt.all = mmp.genes.cnt.all[mmp.use.DB.namelist.for.unip], unip.list = mmp.use.mf.unip.genelist,
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
	mmp.DBgenes.unip.countdf <- bind_rows(mmp.DBgenes.unip.countlist)
	mmp.DBgenes.unip.countdf$DBname <- factor(mmp.DBgenes.unip.countdf$DBname, levels = mmp.use.DB.namelist.for.unip)
	mmp.for.unip.colour <- c(ktDB.mouse.favor.color, ktIT.favor.color)
	names(mmp.for.unip.colour) <- c(ktDB.mouse.namelist, "InterCellDB")
	# plot
	mmp.for.unip.plot <- ggplot(data = mmp.DBgenes.unip.countdf) + 
		geom_col(aes(x = DBname, y = unip.cnt, group = Unip.Keywords, fill = DBname)) + 
		facet_grid(.~Unip.Keywords) + 
		scale_fill_manual(values = mmp.for.unip.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			legend.position = "none")
	# (width = 7.32, height = 4.73)

# [part 004]-sub3: picked GO terms
	mmp.chk.genes.picked.GOlist <- list(
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
	mmp.chk.genes.picked.GOlist <- mmp.chk.genes.picked.GOlist[6:9]
	#
	# get genes
	mmp.chk.genes.picked.goid <- vapply(mmp.chk.genes.picked.GOlist, function(x) {x[2]}, FUN.VALUE = character(1))
	mmp.picked.go.rel.genelist <- lapply(mmp.chk.genes.picked.goid, 
		function(x) {
			Tool.FindGenesFromGO(x, genes.mouse.ref.db, go.mouse.ref.db, 
				go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
		})
	names(mmp.picked.go.rel.genelist) <- vapply(mmp.chk.genes.picked.GOlist, function(x) {x[1]}, FUN.VALUE = character(1))

	# run GO related genes and DB genes intersection
	mmp.use.DB.namelist.picked.GO <- c(ktDB.mouse.namelist, "InterCellDB")
	mmp.DBgenes.GO.countlist <- lapply(seq_along(mmp.genes.cnt.all[mmp.use.DB.namelist.picked.GO]), 
		inside.genes.cnt.all = mmp.genes.cnt.all[mmp.use.DB.namelist.picked.GO], GO.list = mmp.picked.go.rel.genelist, 
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
	mmp.DBgenes.GO.countdf <- bind_rows(mmp.DBgenes.GO.countlist)
	mmp.DBgenes.GO.countdf$DBname <- factor(mmp.DBgenes.GO.countdf$DBname, levels = mmp.use.DB.namelist.picked.GO)
	mmp.DBgenes.GO.countdf$GOterm <- factor(mmp.DBgenes.GO.countdf$GOterm, levels = names(mmp.picked.go.rel.genelist))
	mmp.for.GO.colour <- c(ktDB.mouse.favor.color, ktIT.favor.color)
	names(mmp.for.GO.colour) <- c(ktDB.mouse.namelist, "InterCellDB")
	# plot
	mmp.go.explore.plot <- ggplot(data = mmp.DBgenes.GO.countdf) + 
		geom_col(aes(x = DBname, y = go.cnt, group = GOterm, fill = DBname)) + 
		facet_grid(.~GOterm) + 
		scale_fill_manual(values = mmp.for.GO.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			legend.position = "none")
	# width = 6.57, height = 4.73, follow unip.plot settings






