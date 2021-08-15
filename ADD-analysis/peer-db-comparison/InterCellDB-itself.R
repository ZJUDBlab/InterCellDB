
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database

### --- 1st process ---
## Human part
##  intercelldbDB has gene reference database, but count genes participating gene pairs
	intercelldb.human.genes.res <- unique(c(pairs.human.db$inter.GeneName.A, pairs.human.db$inter.GeneName.B))
	#write.table(intercelldb.human.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-human-intercelldb.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	# Exp sub-DB
	tmp.exp.db <- pairs.human.db[which(pairs.human.db$inter.Experiments.Score != 0), ]
	intercelldb.exp.human.genes.res <- unique(c(tmp.exp.db$inter.GeneName.A, tmp.exp.db$inter.GeneName.B))
	rm(tmp.exp.db)
	# Know sub-DB
	tmp.know.db <- pairs.human.db[intersect(which(pairs.human.db$inter.Experiments.Score == 0), 
		which(pairs.human.db$inter.Database.Score != 0)), ]
	intercelldb.know.human.genes.res <- unique(c(tmp.know.db$inter.GeneName.A, tmp.know.db$inter.GeneName.B))
	rm(tmp.know.db)
	# Pred sub-DB
	tmp.pred.db <- pairs.human.db[intersect(intersect(which(pairs.human.db$inter.Experiments.Score == 0), 
		which(pairs.human.db$inter.Database.Score == 0)), which(pairs.human.db$inter.Predicted.Score != 0)), ]
	intercelldb.pred.human.genes.res <- unique(c(tmp.pred.db$inter.GeneName.A, tmp.pred.db$inter.GeneName.B))
	rm(tmp.pred.db)
	# LR sub-DB
	tmp.pairs.ref <- pairs.human.db#[union(which(pairs.human.db$inter.Experiments.Score != 0), which(pairs.human.db$inter.Database.Score > 700)), ] 
	tmp.pairs.human.lrdb <- FastAlignPairs(tmp.pairs.ref[union(intersect(which(tmp.pairs.ref$inter.GeneName.A %in% sel.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.B %in% sel.ligand.intercelldb.res)),
		intersect(which(tmp.pairs.ref$inter.GeneName.B %in% sel.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.A %in% sel.ligand.intercelldb.res))), ], 4)
	intercelldb.LR.human.genes.res <- unique(c(tmp.pairs.human.lrdb$inter.GeneName.A, tmp.pairs.human.lrdb$inter.GeneName.B))
	rm(tmp.pairs.ref, tmp.pairs.human.lrdb)

## Mouse Part
	intercelldb.mouse.genes.res <- unique(c(pairs.mouse.db$inter.GeneName.A, pairs.mouse.db$inter.GeneName.B))
	#write.table(intercelldb.mouse.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-mouse-intercelldb.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	# Exp sub-DB
	tmp.exp.db <- pairs.mouse.db[which(pairs.mouse.db$inter.Experiments.Score != 0), ]
	intercelldb.exp.mouse.genes.res <- unique(c(tmp.exp.db$inter.GeneName.A, tmp.exp.db$inter.GeneName.B))
	rm(tmp.exp.db)
	# Know sub-DB
	tmp.know.db <- pairs.mouse.db[intersect(which(pairs.mouse.db$inter.Experiments.Score == 0), 
		which(pairs.mouse.db$inter.Database.Score != 0)), ]
	intercelldb.know.mouse.genes.res <- unique(c(tmp.know.db$inter.GeneName.A, tmp.know.db$inter.GeneName.B))
	rm(tmp.know.db)
	# Pred sub-DB
	tmp.pred.db <- pairs.mouse.db[intersect(intersect(which(pairs.mouse.db$inter.Experiments.Score == 0), 
		which(pairs.mouse.db$inter.Database.Score == 0)), which(pairs.mouse.db$inter.Predicted.Score != 0)), ]
	intercelldb.pred.mouse.genes.res <- unique(c(tmp.pred.db$inter.GeneName.A, tmp.pred.db$inter.GeneName.B))
	rm(tmp.pred.db)
	# LR sub-DB
	tmp.pairs.ref <- pairs.mouse.db#[union(which(pairs.mouse.db$inter.Experiments.Score != 0), which(pairs.mouse.db$inter.Database.Score > 700)), ] 
	tmp.pairs.mouse.lrdb <- FastAlignPairs(tmp.pairs.ref[union(intersect(which(tmp.pairs.ref$inter.GeneName.A %in% sel.mouse.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.B %in% sel.mouse.ligand.intercelldb.res)),
		intersect(which(tmp.pairs.ref$inter.GeneName.B %in% sel.mouse.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.A %in% sel.mouse.ligand.intercelldb.res))), ], 4)
	intercelldb.LR.mouse.genes.res <- unique(c(tmp.pairs.mouse.lrdb$inter.GeneName.A, tmp.pairs.mouse.lrdb$inter.GeneName.B))
	rm(tmp.pairs.ref, tmp.pairs.mouse.lrdb)

### --- 2nd process ---
# by function definition

## Human part
# ------
	# option 1.1: get all Cytokine genes, by Uniprot definition on Cytokine
	tmp.cytokine.intercelldb.res <- unique(anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name == "Cytokine"), "Gene.name"])
	write.table(tmp.cytokine.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-intercelldb-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.2: get all Growth factor genes, by Uniprot definition on Growth factor
	tmp.growfactor.intercelldb.res <- unique(anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name == "Growth Factor"), "Gene.name"])
	write.table(tmp.growfactor.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-intercelldb-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.3: get all Receptor genes
	if (FALSE) {  # old definiton, by Uniprot definition on Receptor
		tmp.receptor.keys <- c("Integrin", "G-protein Coupled Receptor", "Host Cell Receptor For Virus Entry", 
			"IgA-binding Protein", "IgE-binding Protein", "IgG-binding Protein", 
			"Photoreceptor Protein", "Receptor", "Retinal Protein")
		dp.tmp.receptor.intercelldb.res <- unique(anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name %in% tmp.receptor.keys), "Gene.name"])
	}
	tmp.receptor.subcellular.score <- c(4, 5)
	sel.receptor.intercelldb.res <- unique(anno.location.human.ref.db[intersect(which(anno.location.human.ref.db$GO.Term.target == "Plasma Membrane"), which(anno.location.human.ref.db$score %in% tmp.receptor.subcellular.score)), "Gene.name"])
	#tmp.receptor.intercelldb.res <- union(tmp.receptor.intercelldb.res, anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name == "Receptor"), "Gene.name"])
	#write.table(tmp.receptor.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/receptor-intercelldb-IT.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.4: get all Ligand genes, by subcellular location on Extracellular Region
	tmp.ligand.subcellular.score <- c(4, 5)
	sel.ligand.intercelldb.res <- unique(anno.location.human.ref.db[intersect(which(anno.location.human.ref.db$GO.Term.target == "Extracellular Region"), which(anno.location.human.ref.db$score %in% tmp.ligand.subcellular.score)), "Gene.name"])
	#write.table(tmp.ligand.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/ligand-intercelldb-IT.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.5: get all Integrin genes, by Uniprot definition on Integrin
	tmp.integrin.intercelldb.res <- unique(anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name == "Integrin"), "Gene.name"])
	write.table(tmp.integrin.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/integrin-intercelldb-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.6: get all G-protein Coupled Receptor, by Uniprot definition on it
	tmp.gpcR.intercelldb.res <- unique(anno.type.human.ref.db[which(anno.type.human.ref.db$Keyword.Name == "G-protein Coupled Receptor"), "Gene.name"])
	write.table(tmp.gpcR.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-intercelldb-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)
# end Human Gene

## Mouse part
# ------
	# option 1.1: get all Cytokine genes, by Uniprot definition on Cytokine
	tmp.mouse.cytokine.intercelldb.res <- unique(anno.type.mouse.ref.db[which(anno.type.mouse.ref.db$Keyword.Name == "Cytokine"), "Gene.name"])
	write.table(tmp.mouse.cytokine.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/cytokine-intercelldb-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.2: get all Growth factor genes, by Uniprot definition on Growth factor
	tmp.mouse.growfactor.intercelldb.res <- unique(anno.type.mouse.ref.db[which(anno.type.mouse.ref.db$Keyword.Name == "Growth Factor"), "Gene.name"])
	write.table(tmp.mouse.growfactor.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/growfactor-intercelldb-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.3: get all Receptor genes
	tmp.mouse.receptor.subcellular.score <- c(4,5)
	sel.mouse.receptor.intercelldb.res <- unique(anno.location.mouse.ref.db[intersect(which(anno.location.mouse.ref.db$GO.Term.target == "Plasma Membrane"), which(anno.location.mouse.ref.db$score %in% tmp.mouse.receptor.subcellular.score)), "Gene.name"])
	#sel.mouse.receptor.intercelldb.res <- union(tmp.mouse.receptor.intercelldb.res, anno.type.mouse.ref.db[which(anno.type.mouse.ref.db$Keyword.Name == "Receptor"), "Gene.name"])
	#write.table(tmp.mouse.receptor.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/receptor-intercelldb-mouse.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.4: get all Ligand genes, by subcellular location on Extracellular Region
	tmp.mouse.ligand.subcellular.score <- c(4,5)
	sel.mouse.ligand.intercelldb.res <- unique(anno.location.mouse.ref.db[intersect(which(anno.location.mouse.ref.db$GO.Term.target == "Extracellular Region"), which(anno.location.mouse.ref.db$score %in% tmp.mouse.ligand.subcellular.score)), "Gene.name"])
	#write.table(tmp.mouse.ligand.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/ligand-intercelldb-mouse.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.5: get all Integrin genes, by Uniprot definition on Integrin
	tmp.mouse.integrin.intercelldb.res <- unique(anno.type.mouse.ref.db[which(anno.type.mouse.ref.db$Keyword.Name == "Integrin"), "Gene.name"])
	write.table(tmp.mouse.integrin.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/integrin-intercelldb-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.6: get all G-protein Coupled Receptor, by Uniprot definition on it
	tmp.mouse.gpcR.intercelldb.res <- unique(anno.type.mouse.ref.db[which(anno.type.mouse.ref.db$Keyword.Name == "G-protein Coupled Receptor"), "Gene.name"])
	write.table(tmp.mouse.gpcR.intercelldb.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/gpcR-intercelldb-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)
# end Mouse Gene


### --- 3rd process ---
# by definition of subcellular locations and molecular functions
# - subcellular locations -
splits.intercelldb.scloc <- levels(factor(anno.location.human.ref.db$GO.Term.target))
# Human get 13 kinds of scloc

# - molecular function -
splits.intercelldb.mfunc <- levels(factor(anno.type.human.ref.db$Keyword.Name))
# Human get 131 kinds of mfunc





##### For Describing PAIRS
### process <0> defining InterCellDB sub-DBs

## Human part
	# Experiments
	pairs.human.experiments.ref <- pairs.human.db[which(pairs.human.db$inter.Experiments.Score != 0), ] 
	pairs.human.experiments.ref <- FastAlignPairs(pairs.human.experiments.ref, 4)
	# position-align
	pos.intercelldb.pairs.human.experiments <- list(gene.A = pairs.human.experiments.ref$inter.GeneName.A, gene.B = pairs.human.experiments.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.human.experiments <- paste(pairs.human.experiments.ref$inter.GeneName.A, pairs.human.experiments.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.human.experiments.ref)
	#write.table(ref.intercelldb.pairs.human.experiments, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ref-pairs-intercelldb-human-experiments.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	# Experiments -sub
		# high
		pairs.human.experiments.highconf <- pairs.human.db[which(pairs.human.db$inter.Experiments.Score > 700), ] 
		pairs.human.experiments.highconf <- FastAlignPairs(pairs.human.experiments.highconf, 4)
		highconf.intercelldb.pairs.human.experiments <- paste(pairs.human.experiments.highconf$inter.GeneName.A, pairs.human.experiments.highconf$inter.GeneName.B, sep = kgp.cut)
		rm(pairs.human.experiments.highconf)
		# middle
		pairs.human.experiments.midconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Experiments.Score <= 700),
			which(pairs.human.db$inter.Experiments.Score > 400)), ] 
		pairs.human.experiments.midconf <- FastAlignPairs(pairs.human.experiments.midconf, 4)
		midconf.intercelldb.pairs.human.experiments <- paste(pairs.human.experiments.midconf$inter.GeneName.A, pairs.human.experiments.midconf$inter.GeneName.B, sep = kgp.cut)
		rm(pairs.human.experiments.midconf)
		# low
		pairs.human.experiments.lowconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Experiments.Score <= 400),
			which(pairs.human.db$inter.Experiments.Score > 0)), ]
		pairs.human.experiments.lowconf <- FastAlignPairs(pairs.human.experiments.lowconf, 4)
		lowconf.intercelldb.pairs.human.experiments <- paste(pairs.human.experiments.lowconf$inter.GeneName.A, pairs.human.experiments.lowconf$inter.GeneName.B, sep = kgp.cut)
		rm(pairs.human.experiments.lowconf)


	# from Database
	pairs.human.knowledge.ref <- pairs.human.db[intersect(which(pairs.human.db$inter.Experiments.Score == 0), 
		which(pairs.human.db$inter.Database.Score != 0)), ] 
	pairs.human.knowledge.ref <- FastAlignPairs(pairs.human.knowledge.ref, 4)
	# position-align
	pos.intercelldb.pairs.human.knowledge <- list(gene.A = pairs.human.knowledge.ref$inter.GeneName.A, gene.B = pairs.human.knowledge.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.human.knowledge <- paste(pairs.human.knowledge.ref$inter.GeneName.A, pairs.human.knowledge.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.human.knowledge.ref)
	#write.table(ref.intercelldb.pairs.human.knowledge, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ref-pairs-intercelldb-human-knowledge.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	# Database -sub
		# high
		pairs.human.knowledge.highconf <- pairs.human.db[which(pairs.human.db$inter.Database.Score > 700), ] 
		pairs.human.knowledge.highconf <- FastAlignPairs(pairs.human.knowledge.highconf, 4)
		highconf.intercelldb.pairs.human.knowledge <- paste(pairs.human.knowledge.highconf$inter.GeneName.A, pairs.human.knowledge.highconf$inter.GeneName.B, sep = kgp.cut)
		highconf.intercelldb.pairs.human.knowledge <- intersect(highconf.intercelldb.pairs.human.knowledge, ref.intercelldb.pairs.human.knowledge)
		rm(pairs.human.knowledge.highconf)
		# high
		pairs.human.knowledge.midconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Database.Score <= 700), 
			which(pairs.human.db$inter.Database.Score > 400)), ] 
		pairs.human.knowledge.midconf <- FastAlignPairs(pairs.human.knowledge.midconf, 4)
		midconf.intercelldb.pairs.human.knowledge <- paste(pairs.human.knowledge.midconf$inter.GeneName.A, pairs.human.knowledge.midconf$inter.GeneName.B, sep = kgp.cut)
		midconf.intercelldb.pairs.human.knowledge <- intersect(midconf.intercelldb.pairs.human.knowledge, ref.intercelldb.pairs.human.knowledge)
		rm(pairs.human.knowledge.midconf)	
		# high
		pairs.human.knowledge.lowconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Database.Score <= 400), 
			which(pairs.human.db$inter.Database.Score > 0)), ] 
		pairs.human.knowledge.lowconf <- FastAlignPairs(pairs.human.knowledge.lowconf, 4)
		lowconf.intercelldb.pairs.human.knowledge <- paste(pairs.human.knowledge.lowconf$inter.GeneName.A, pairs.human.knowledge.lowconf$inter.GeneName.B, sep = kgp.cut)
		lowconf.intercelldb.pairs.human.knowledge <- intersect(lowconf.intercelldb.pairs.human.knowledge, ref.intercelldb.pairs.human.knowledge)
		rm(pairs.human.knowledge.lowconf)


	# by Prediction
	pairs.human.prediction.ref <- pairs.human.db[intersect(intersect(which(pairs.human.db$inter.Experiments.Score == 0), 
		which(pairs.human.db$inter.Database.Score == 0)), which(pairs.human.db$inter.Predicted.Score != 0)), ]
	pairs.human.prediction.ref <- FastAlignPairs(pairs.human.prediction.ref, 4)
	# position-align
	pos.intercelldb.pairs.human.prediction <- list(gene.A = pairs.human.prediction.ref$inter.GeneName.A, gene.B = pairs.human.prediction.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.human.prediction <- paste(pairs.human.prediction.ref$inter.GeneName.A, pairs.human.prediction.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.human.prediction.ref)
	#write.table(ref.intercelldb.pairs.human.prediction, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ref-pairs-intercelldb-human-prediction.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	# Prediction -sub
		# high
		pairs.human.prediction.highconf <- pairs.human.db[which(pairs.human.db$inter.Predicted.Score > 700), ] 
		pairs.human.prediction.highconf <- FastAlignPairs(pairs.human.prediction.highconf, 4)
		highconf.intercelldb.pairs.human.prediction <- paste(pairs.human.prediction.highconf$inter.GeneName.A, pairs.human.prediction.highconf$inter.GeneName.B, sep = kgp.cut)
		highconf.intercelldb.pairs.human.prediction <- intersect(highconf.intercelldb.pairs.human.prediction, ref.intercelldb.pairs.human.prediction)
		rm(pairs.human.prediction.highconf)
		# high
		pairs.human.prediction.midconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Predicted.Score <= 700), 
			which(pairs.human.db$inter.Predicted.Score > 400)), ] 
		pairs.human.prediction.midconf <- FastAlignPairs(pairs.human.prediction.midconf, 4)
		midconf.intercelldb.pairs.human.prediction <- paste(pairs.human.prediction.midconf$inter.GeneName.A, pairs.human.prediction.midconf$inter.GeneName.B, sep = kgp.cut)
		midconf.intercelldb.pairs.human.prediction <- intersect(midconf.intercelldb.pairs.human.prediction, ref.intercelldb.pairs.human.prediction)
		rm(pairs.human.prediction.midconf)	
		# high
		pairs.human.prediction.lowconf <- pairs.human.db[intersect(which(pairs.human.db$inter.Predicted.Score <= 400), 
			which(pairs.human.db$inter.Predicted.Score > 0)), ] 
		pairs.human.prediction.lowconf <- FastAlignPairs(pairs.human.prediction.lowconf, 4)
		lowconf.intercelldb.pairs.human.prediction <- paste(pairs.human.prediction.lowconf$inter.GeneName.A, pairs.human.prediction.lowconf$inter.GeneName.B, sep = kgp.cut)
		lowconf.intercelldb.pairs.human.prediction <- intersect(lowconf.intercelldb.pairs.human.prediction, ref.intercelldb.pairs.human.prediction)
		rm(pairs.human.prediction.lowconf)

	# LR
	tmp.pairs.ref <- pairs.human.db#[union(which(pairs.human.db$inter.Experiments.Score != 0), which(pairs.human.db$inter.Database.Score > 700)), ] 
	tmp.pairs.human.lrdb <- FastAlignPairs(tmp.pairs.ref[union(intersect(which(tmp.pairs.ref$inter.GeneName.A %in% sel.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.B %in% sel.ligand.intercelldb.res)),
		intersect(which(tmp.pairs.ref$inter.GeneName.B %in% sel.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.A %in% sel.ligand.intercelldb.res))), ], 4)
	# position-align
	pos.intercelldb.pairs.human.LR <- list(gene.A = tmp.pairs.human.lrdb$inter.GeneName.A, gene.B = tmp.pairs.human.lrdb$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.human.LR <- paste(tmp.pairs.human.lrdb$inter.GeneName.A, tmp.pairs.human.lrdb$inter.GeneName.B, sep = kgp.cut)
	rm(tmp.pairs.human.lrdb, tmp.pairs.ref)


	# ALL
	pairs.human.ref <- FastAlignPairs(pairs.human.db, 4)
	# position-align
	pos.intercelldb.pairs.human <- list(gene.A = pairs.human.ref$inter.GeneName.A, gene.B = pairs.human.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.human <- paste(pairs.human.ref$inter.GeneName.A, pairs.human.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.human.ref)
	#write.table(ref.intercelldb.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ref-pairs-intercelldb-human-all.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# Action mode
	ref.intercelldb.act.mode.pairs.human <- lapply(ktDB.act.mode, ref.action.db = actions.human.ref.db,
		function(x, ref.action.db) {
			tmp.ret.db <- unique(FastAlignPairs(ref.action.db[which(ref.action.db$mode == x), 1:4], 4))  # only ID and Gene.name
			paste(tmp.ret.db$inter.GeneName.A, tmp.ret.db$inter.GeneName.B, sep = kgp.cut)
		})
	names(ref.intercelldb.act.mode.pairs.human) <- ktDB.act.mode

	# Action effect
	tmp.is.directional <- which(actions.human.ref.db$is_directional == 't')
	tmp.not.directional <- which(actions.human.ref.db$is_directional == 'f')
	tmp.is.aisacting <- which(actions.human.ref.db$a_is_acting == 't')
	tmp.not.aisacting <- which(actions.human.ref.db$a_is_acting == 'f')
	ref.intercelldb.act.effect.pairs.human <- lapply(ktDB.act.effect, 
		ref.action.db = actions.human.ref.db, 
		is.directional = tmp.is.directional, not.directional = tmp.not.directional,
		is.aisacting = tmp.is.aisacting, not.aisacting = tmp.not.aisacting, 
		function(x, ref.action.db, is.directional, not.directional, is.aisacting, not.aisacting) {
			tmp.inds <- switch(x,
				"positive" = {
					Reduce(intersect, list(is.directional, is.aisacting, 
						which(ref.action.db$action == "activation")))
				},
				"negative" = {
					Reduce(intersect, list(is.directional, is.aisacting, 
						which(ref.action.db$action == "inhibition")))
				},
				"unspecified" = {
					Reduce(intersect, list(is.directional, is.aisacting, 
						which(ref.action.db$action == "")))
				},
				"undirected" = {
					intersect(not.directional, not.aisacting)
				},
				stop("Undefined action effects!"))
			tmp.ret.db <- unique(FastAlignPairs(ref.action.db[tmp.inds, 1:4], 4))  # only ID and Gene.name
			paste(tmp.ret.db$inter.GeneName.A, tmp.ret.db$inter.GeneName.B, sep = kgp.cut)
		})
	names(ref.intercelldb.act.effect.pairs.human) <- ktDB.act.effect

# --- end ---




## Mouse part
	# Experiments
	pairs.mouse.experiments.ref <- pairs.mouse.db[which(pairs.mouse.db$inter.Experiments.Score != 0), ] 
	pairs.mouse.experiments.ref <- FastAlignPairs(pairs.mouse.experiments.ref, 4)
	# position-align
	pos.intercelldb.pairs.mouse.experiments <- list(gene.A = pairs.mouse.experiments.ref$inter.GeneName.A, gene.B = pairs.mouse.experiments.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.mouse.experiments <- paste(pairs.mouse.experiments.ref$inter.GeneName.A, pairs.mouse.experiments.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.mouse.experiments.ref)
	#write.table(ref.intercelldb.pairs.mouse.experiments, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/ref-pairs-intercelldb-mouse-experiments.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# from Database
	pairs.mouse.knowledge.ref <- pairs.mouse.db[intersect(which(pairs.mouse.db$inter.Experiments.Score == 0), 
		which(pairs.mouse.db$inter.Database.Score != 0)), ] 
	pairs.mouse.knowledge.ref <- FastAlignPairs(pairs.mouse.knowledge.ref, 4)
	# position-align
	pos.intercelldb.pairs.mouse.knowledge <- list(gene.A = pairs.mouse.knowledge.ref$inter.GeneName.A, gene.B = pairs.mouse.knowledge.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.mouse.knowledge <- paste(pairs.mouse.knowledge.ref$inter.GeneName.A, pairs.mouse.knowledge.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.mouse.knowledge.ref)
	#write.table(ref.intercelldb.pairs.mouse.knowledge, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/ref-pairs-intercelldb-mouse-knowledge.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# by Prediction
	pairs.mouse.prediction.ref <- pairs.mouse.db[intersect(intersect(which(pairs.mouse.db$inter.Experiments.Score == 0), 
		which(pairs.mouse.db$inter.Database.Score == 0)), which(pairs.mouse.db$inter.Predicted.Score != 0)), ]
	pairs.mouse.prediction.ref <- FastAlignPairs(pairs.mouse.prediction.ref, 4)
	# position-align
	pos.intercelldb.pairs.mouse.prediction <- list(gene.A = pairs.mouse.prediction.ref$inter.GeneName.A, gene.B = pairs.mouse.prediction.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.mouse.prediction <- paste(pairs.mouse.prediction.ref$inter.GeneName.A, pairs.mouse.prediction.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.mouse.prediction.ref)
	#write.table(ref.intercelldb.pairs.mouse.prediction, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/ref-pairs-intercelldb-mouse-prediction.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# LR
	tmp.pairs.ref <- pairs.mouse.db#[union(which(pairs.mouse.db$inter.Experiments.Score != 0), which(pairs.mouse.db$inter.Database.Score > 700)), ] 
	tmp.pairs.mouse.lrdb <- FastAlignPairs(tmp.pairs.ref[union(intersect(which(tmp.pairs.ref$inter.GeneName.A %in% sel.mouse.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.B %in% sel.mouse.ligand.intercelldb.res)),
		intersect(which(tmp.pairs.ref$inter.GeneName.B %in% sel.mouse.receptor.intercelldb.res), which(tmp.pairs.ref$inter.GeneName.A %in% sel.mouse.ligand.intercelldb.res))), ], 4)
	# position-align
	pos.intercelldb.pairs.mouse.LR <- list(gene.A = tmp.pairs.mouse.lrdb$inter.GeneName.A, gene.B = tmp.pairs.mouse.lrdb$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.mouse.LR <- paste(tmp.pairs.mouse.lrdb$inter.GeneName.A, tmp.pairs.mouse.lrdb$inter.GeneName.B, sep = kgp.cut)
	rm(tmp.pairs.mouse.lrdb, tmp.pairs.ref)

	# ALL
	pairs.mouse.ref <- FastAlignPairs(pairs.mouse.db, 4)
	# position-align
	pos.intercelldb.pairs.mouse <- list(gene.A = pairs.mouse.ref$inter.GeneName.A, gene.B = pairs.mouse.ref$inter.GeneName.B)
	# pasted style
	ref.intercelldb.pairs.mouse <- paste(pairs.mouse.ref$inter.GeneName.A, pairs.mouse.ref$inter.GeneName.B, sep = kgp.cut)
	rm(pairs.mouse.ref)
	#write.table(ref.intercelldb.pairs.mouse, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/ref-pairs-intercelldb-mouse-all.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

