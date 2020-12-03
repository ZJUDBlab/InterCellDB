

# get other pairs
other.pairs.names.C <- setdiff(grep(paste0("^", this.cluster.C), interact.pairs.acted$name.allpairs, value = TRUE), involved.clusters.pair)
other.pairs.names.D <- setdiff(grep(paste0(this.cluster.D, "$"), interact.pairs.acted$name.allpairs, value = TRUE), involved.clusters.pair)


FindSpecialGenesInOnepairCluster <- function(
	interact.pairs.acted,
	target.gene.pairs.df,  # get from VEinfos [TODO]
	involved.clusters.pair, # = paste("cluster.C", "cluster.D", sep = kClustersSplit),
	to.cmp.clusters.pairs = character(),
	to.cmp.clusters.pairs.sel.strategy = "inter-cross",
	sep.inside.gene.pair = "->",
	...
) {
	this.gp.sep <- sep.inside.gene.pair
	# calculate the multiply of fold change of interacting genes
	inside.c.short.interacts <- function(
		tmp.pairs,
		tmp.sep
	) {
		paste(tmp.pairs[, "inter.GeneName.A"], tmp.pairs[, "inter.GeneName.B"], sep = tmp.sep)
	}

	# all pairs
	all.pairs.names <- interact.pairs.acted$name.allpairs
	all.pairs.interacts <- interact.pairs.acted$data.allpairs
	# this pair
	this.pair.name <- involved.clusters.pair
	this.clusters.separate <- strsplit(this.pair.name, split = kClustersSplit, fixed = TRUE)[[1]]
	this.pair.interacts <- inside.c.short.interacts(target.gene.pairs.df, this.gp.sep)

	# to compare pairs
	if (length(to.cmp.clusters.pairs) == 0) {  # only if no list is given then use the pre-defined strategy
		if (to.cmp.clusters.pairs.sel.strategy == "inter-cross") {
			other.pairs.names.C <- setdiff(grep(paste0("^", this.clusters.separate[1]), all.pairs.names, value = TRUE), this.pair.name)
			other.pairs.names.D <- setdiff(grep(paste0(this.clusters.separate[2], "$"), all.pairs.names, value = TRUE), this.pair.name)
			to.cmp.clusters.pairs <- unique(c(other.pairs.names.C, other.pairs.names.D))
		} else {
			if (to.cmp.clusters.pairs.sel.strategy == "all-other") {
				to.cmp.clusters.pairs <- setdiff(all.pairs.names, this.pair.name)
			} else {
				stop("Undefined pre-defined strategy used: ", to.cmp.clusters.pairs.sel.strategy)
			}
		}
	}  # else use the directly specified clusters

	# generate gene pairs to compare
	other.packed.gene.pairs.df <- lapply(to.cmp.clusters.pairs, all.pairs.interacts = all.pairs.interacts, 
		function(x, all.pairs.interacts) { all.pairs.interacts[[x]] })
	names(other.packed.gene.pairs.df) <- to.cmp.clusters.pairs
	# get the short interacts
	other.packed.pair.interacts <- lapply(other.packed.gene.pairs.df, sep = this.gp.sep, 
		function(x, sep) { inside.c.short.interacts(x, this.gp.sep) })
	

	# get diff pairs

	

}

