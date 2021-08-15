
# This file is to explore the hierarchy tree of Uniprot Keywords 
# This file works with InterCellDB::anno.type.***.db

library(InterCellDB)
# for construct class
library(methods)
# for construct hierarchy tree
library(ggraph)
library(igraph)
library(tidyverse)

# fetch KeyWord list (mouse and human are of little differenct)
tmp.human.types <- levels(factor(anno.type.human.ref.db$Keyword.Name))
tmp.mouse.types <- levels(factor(anno.type.mouse.ref.db$Keyword.Name))
if (FALSE) {  # check overlaps
	length(which(tmp.mouse.types %in% tmp.human.types))  # 125
	setdiff(tmp.mouse.types, tmp.human.types)  # only 1 keyword, "Pheromone"
}

# As human get more Keywords, use human data to generate the hierarchy tree


### the process
uniprot.inc.all <- unique(union(tmp.human.types, tmp.mouse.types))
uniprot.inc.genelist <- lapply(uniprot.inc.all, ref.type.db = anno.type.human.ref.db,
	function(x, ref.type.db) {
		ref.type.db[(ref.type.db$Keyword.Name == x), "Gene.name"]
		})
## sequentially construct the tree
# set class, each node represents one Keyword
UnipNode <- setClass("UnipNode",
	slots = c(node.it = "character",
		node.parents = "character",
		node.children = "character"))


# least method (N^2) compare to all other and direct get UnipNode
unip.node.list <- lapply(seq_along(uniprot.inc.all),
	Keyref.name = uniprot.inc.all, Keyref.belongs = uniprot.inc.genelist, 
	unstable.offset.percentage = .05,  # details put under this function
	number.limit.overmatching = 3,  # genes under xxx, may get too many parents, use this to cut those relations and leave it to manual work
	function(x, Keyref.name, Keyref.belongs, number.limit.overmatching, unstable.offset.percentage) {
		this.ret <- new(Class = "UnipNode", node.it = Keyref.name[x])
		# find parents
		this.is.parents <- vector(mode = "logical", length = length(Keyref.name))
		for (i in seq_along(Keyref.belongs)) {
			tmp.matching <- length(which(Keyref.belongs[[x]] %in% Keyref.belongs[[i]]))
			tmp.match.interval <- length(Keyref.belongs[[x]]) * (unstable.offset.percentage * c(-1, 1) + 1)
			tmp.match.interval <- c(floor(tmp.match.interval[1]), ceiling(tmp.match.interval[2]))
			this.is.parents[i] <- (length(Keyref.belongs[[x]]) != 0) && ((tmp.matching >= tmp.match.interval[1]) && (tmp.matching <= tmp.match.interval[2]))
		}



		# [TODO] number.limit.overmatching



		this.ret@node.parents <- setdiff(Keyref.name[this.is.parents], Keyref.name[x])
		# find children
		this.is.children <- vector(mode = "logical", length = length(Keyref.name))
		for (i in seq_along(Keyref.belongs)) {
			tmp.matching <- length(which(Keyref.belongs[[i]] %in% Keyref.belongs[[x]]))
			tmp.match.interval <- length(Keyref.belongs[[i]]) * (unstable.offset.percentage * c(-1, 1) + 1)
			tmp.match.interval <- c(floor(tmp.match.interval[1]), ceiling(tmp.match.interval[2]))
			this.is.children[i] <- (length(Keyref.belongs[[x]]) != 0) && ((tmp.matching >= tmp.match.interval[1]) && (tmp.matching <= tmp.match.interval[2]))
		}








		this.ret@node.children <- setdiff(Keyref.name[this.is.children], Keyref.name[x])
		return(this.ret)
	})
# unstable.offset Explanation
# In data, see "IgE-binding Protein" %in% "Receptor" only get 1 gene unmatched, though total 5 genes
# however, it is biologically belonging relationship.
# So this parameter is invented to give 5% unmatching part to make nicer automatical hierarchy construction.


# generate from-to format for drawing hierarchy plot
unip.edges <- lapply(unip.node.list, 
	function(x) {
		this.parents.to.it <- data.frame(from = x@node.parents, to = rep(x@node.it, times = length(x@node.parents)), stringsAsFactors = FALSE)
		this.it.to.children <- data.frame(from = rep(x@node.it, times = length(x@node.children)), to = x@node.children, stringsAsFactors = FALSE)
		# no parents nodes, will be applied with dummy parents
		this.dummy <- this.parents.to.it[integer(0), ]
		if (nrow(this.parents.to.it) == 0) {
			this.dummy <- data.frame(from = "dummy.root", to = x@node.it, stringsAsFactors = FALSE)
		}
		Reduce(rbind, list(this.parents.to.it, this.it.to.children, this.dummy))
	})
unip.edges <- dplyr::bind_rows(unip.edges)
unip.edges <- unique(unip.edges)



## plot the hierarchy 
# Create a graph object 
mygraph <- graph_from_data_frame(unip.edges)
V(mygraph)$node_label <- attr(V(mygraph), "names")
 
# Basic tree
ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_bend(colour = "grey", alpha = 0.8) +
  geom_node_point(colour = "red") +
  geom_node_text(aes(label = node_label), size = 3, angle = 90) + 
  #scale_x_continuous(expand = expansion(mult = 0.5)) + 
  scale_y_continuous(expand = expansion(mult = 0.5)) + 
  theme_void()

## Final Collection
# manual collect all relations by this function
# and compare with functional merging by Xiaotao
# ---
# result generated





