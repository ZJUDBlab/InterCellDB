actions.to.op.ref.db <- actions.human.ref.db
# unique on gene pairs
ref.slim.db <- unique(FastAlignPairs(actions.to.op.ref.db[, 1:4], 4))
ref.slim.align.conv <- paste(ref.slim.db$inter.GeneName.A, ref.slim.db$inter.GeneName.B, sep = ">")
ref.slim.align.rev <- paste(ref.slim.db$inter.GeneName.B, ref.slim.db$inter.GeneName.A, sep = ">")
# original database
ref.pairs.align <- paste(actions.to.op.ref.db$inter.GeneName.A, actions.to.op.ref.db$inter.GeneName.B,
	sep = ">")
# 
ref.pairs.inds.collect <- tapply(seq_along(ref.pairs.align), ref.pairs.align, function(x) {x})

library(future.apply)
# new way to get unique actions
#prog.bar.c1 <- progress::progress_bar$new(total = length(ref.slim.align.conv))
#prog.bar.c1$tick(0)
inds.actions.subset <- future_lapply(seq_along(ref.slim.align.conv), 
	slim.conv = ref.slim.align.conv, slim.rev = ref.slim.align.rev,
	ref.inds.collect = ref.pairs.inds.collect, ref.actions.db = actions.to.op.ref.db, 
	function(x, slim.conv, slim.rev, ref.inds.collect, ref.actions.db) {
		one.gp.actions.db <- ref.actions.db[c(ref.inds.collect[[slim.conv[x]]], ref.inds.collect[[slim.rev[x]]]), ]
		inds.one.gp.conv <- seq_along(ref.inds.collect[[slim.conv[x]]])
		inds.one.gp.rev <- seq_along(ref.inds.collect[[slim.rev[x]]]) + length(ref.inds.collect[[slim.conv[x]]])
		# remove all is_directional = t but a_is_acting = f ones. Those are redudant informations
		#inds.tf <- intersect(which(one.gp.actions.db$is_directional == 't'), which(one.gp.actions.db$a_is_acting == 'f'))
		#one.gp.actions.db <- one.gp.actions.db[setdiff(seq_len(nrow(one.gp.actions.db)), inds.tf), ]
		one.gp.gene.ab <- strsplit(slim.conv[x], split = ">", fixed = TRUE)[[1]]
		one.gp.mode.act <- paste(one.gp.actions.db$mode, one.gp.actions.db$action, sep = "~")
		one.gp.mode.act.sels <- tapply(seq_along(one.gp.mode.act), one.gp.mode.act, function(x) {x})
		tmp.res.mode.act <- lapply(one.gp.mode.act.sels, one.gp.ref.db = one.gp.actions.db, 
			one.gp.gene.ab = one.gp.gene.ab, inds.one.gp.conv = inds.one.gp.conv, inds.one.gp.rev = inds.one.gp.rev, 
			function(x, one.gp.ref.db, one.gp.gene.ab, inds.one.gp.conv, inds.one.gp.rev) {
				ret.orig.inds <- integer()
				# conv part
				tmp.inds.conv <- intersect(x, inds.one.gp.conv)
				if (length(tmp.inds.conv) > 0) {
					tmp.conv.db <- one.gp.ref.db[tmp.inds.conv, ]
					# with t-f, f-t conditions removed
					inds.topsel <- intersect(which(tmp.conv.db$is_directional == 't'), which(tmp.conv.db$a_is_acting == 't'))
					inds.lastsel <- intersect(which(tmp.conv.db$is_directional == 'f'), which(tmp.conv.db$a_is_acting == 'f'))
					if (length(inds.topsel) > 0) {
						inds.sel.ret <- inds.topsel
					} else {  # inds.lastsel must > 0 when all other not exists
							inds.sel.ret <- inds.lastsel
					}
					if (length(inds.sel.ret) > 1) {
						tmp.db <- tmp.conv.db[inds.sel.ret, ]
						ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.db)[which(tmp.db$score == max(tmp.db$score))]))
					} else {
						if (length(inds.sel.ret) > 0) {
							ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.conv.db)[inds.sel.ret]))
						}
					}
				}
				
				# rev part
				tmp.inds.rev <- intersect(x, inds.one.gp.rev)
				if (length(tmp.inds.rev) > 0) {
					tmp.rev.db <- one.gp.ref.db[tmp.inds.rev, ]
					# with t-f, f-t conditions removed
					inds.topsel <- intersect(which(tmp.rev.db$is_directional == 't'), which(tmp.rev.db$a_is_acting == 't'))
					inds.lastsel <- intersect(which(tmp.rev.db$is_directional == 'f'), which(tmp.rev.db$a_is_acting == 'f'))
					if (length(inds.topsel) > 0) {
						inds.sel.ret <- inds.topsel
					} else {  # inds.lastsel must > 0 when all other not exists
							inds.sel.ret <- inds.lastsel
					}
					if (length(inds.sel.ret) > 1) {
						tmp.db <- tmp.rev.db[inds.sel.ret, ]
						ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.db)[which(tmp.db$score == max(tmp.db$score))]))
					} else {
						if (length(inds.sel.ret) > 0) {
							ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.rev.db)[inds.sel.ret]))
						}
					}
				}
				ret.orig.inds
			})
		#prog.bar.c1$tick()
		as.integer(unlist(tmp.res.mode.act))
	})
timestamp()
# fetch result
inds.result.actions.subset <- as.integer(unlist(inds.actions.subset))
actions.uq.mouse.db <- actions.to.op.ref.db[which(rownames(actions.to.op.ref.db) %in% inds.result.actions.subset),]
# done, save RDS



# old way to get unique actions
if (FALSE) {
	prog.bar.c1 <- progress::progress_bar$new(total = length(ref.slim.align.conv))
	prog.bar.c1$tick(0)
	inds.actions.subset <- lapply(seq_along(ref.slim.align.conv), 
	slim.conv = ref.slim.align.conv, slim.rev = ref.slim.align.rev,
	ref.pairs.align = ref.pairs.align, ref.actions.db = actions.human.ref.db, 
	function(x, slim.conv, slim.rev, ref.pairs.align, ref.actions.db) {
		one.gp.actions.db <- ref.actions.db[union(which(ref.pairs.align == slim.conv[x]),
			which(ref.pairs.align == slim.rev[x])), ]
		one.gp.gene.ab <- strsplit(slim.conv[x], split = ">", fixed = TRUE)[[1]]
		one.gp.mode.sels <- unique(one.gp.actions.db$mode)
		one.gp.act.sels <- unique(one.gp.actions.db$action)
		ret.orig.inds <- integer()
		for (i in one.gp.mode.sels) {
			for (j in one.gp.act.sels) {
				for (k in 1:2) {  # k = 1, check conv, k = 2, check rev
					inds.chk.exist <- intersect(which(one.gp.actions.db$mode == i), which(one.gp.actions.db$action == j))
					if (length(inds.chk.exist) == 0) {
						next
					}
					inside.gp.gene.ab <- one.gp.gene.ab
					if (k == 2) {
						inside.gp.gene.ab <- c(inside.gp.gene.ab[2], inside.gp.gene.ab[1])
					}
					inds.it <- intersect(which(one.gp.actions.db$inter.GeneName.A == inside.gp.gene.ab[1]),
							which(one.gp.actions.db$inter.GeneName.B == inside.gp.gene.ab[2]))
					inds.it <- intersect(inds.it, inds.chk.exist)
					one.gp.sel.sub.db <- one.gp.actions.db[inds.it, ]
					if (length(inds.it) == 0) {
						stop("Unexpected Error!")
					} else {
						inds.topsel <- intersect(which(one.gp.sel.sub.db$is_directional == 't'), which(one.gp.sel.sub.db$a_is_acting == 't'))
						inds.secsel <- intersect(which(one.gp.sel.sub.db$is_directional == 't'), which(one.gp.sel.sub.db$a_is_acting == 'f'))
						inds.lastsel <- intersect(which(one.gp.sel.sub.db$is_directional == 'f'), which(one.gp.sel.sub.db$a_is_acting == 'f'))
						inds.error.chk <- intersect(which(one.gp.sel.sub.db$is_directional == 'f'), which(one.gp.sel.sub.db$a_is_acting == 't'))
						if (length(inds.error.chk) > 0) {
							stop("Unexpected database items, with is_directional == f, but a_is_acting == t!")
						}
						inds.sel.ret <- integer(0)
						if (length(inds.topsel) > 0) {
							inds.sel.ret <- inds.topsel
						} else {
							if (length(inds.secsel) > 0) {
								# check if the rev has the topsel for the same mode and action
								tmp.inds.it.rev <- Reduce(intersect, list(which(one.gp.actions.db$inter.GeneName.A == inside.gp.gene.ab[2]), 
									which(one.gp.actions.db$inter.GeneName.B == inside.gp.gene.ab[1]),
									which(one.gp.actions.db$is_directional == 't'),
									which(one.gp.actions.db$a_is_acting == 't')))
								tmp.inds.it.rev <- intersect(tmp.inds.it.rev, inds.chk.exist)  # add that mode and action the same
								if (length(tmp.inds.it.rev) == 0) {  # rev has no explicit definitions, so keep this records
									inds.sel.ret <- inds.secsel
								}
							} else {  # inds.lastsel must > 0 when all other not exists
								inds.sel.ret <- inds.lastsel
							}
						}
						if (length(inds.sel.ret) > 1) {
							message(paste0(slim.conv[x], " in mode: ", i, ", action: ", j, ", has multiple score!"))
							tmp.db <- one.gp.sel.sub.db[inds.sel.ret, ]
							tmp.db <- tmp.db[order(tmp.db$score, decreasing = TRUE), ]
							ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.db[1])))
						} else {
							if (length(inds.sel.ret) > 0) {
								ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(one.gp.sel.sub.db[inds.sel.ret, ])))
							}
						}
					}
				}
			}
		}
		prog.bar.c1$tick()
		ret.orig.inds
	})
}

# modified version to detect t-f pairs
#
#
ref.tf.db <- unique(FastAlignPairs(actions.human.ref.db[Reduce(intersect, list(
	which(actions.human.ref.db$is_directional == 't'),
	which(actions.human.ref.db$a_is_acting == 'f'))), 1:4], 4))
ref.tf.align.conv <- paste(ref.tf.db$inter.GeneName.A, ref.tf.db$inter.GeneName.B, sep = ">")
ref.tf.align.rev <- paste(ref.tf.db$inter.GeneName.B, ref.tf.db$inter.GeneName.A, sep = ">")
# original database
ref.pairs.align <- paste(actions.human.ref.db$inter.GeneName.A, actions.human.ref.db$inter.GeneName.B,
	sep = ">")

inds.to.explore <- sample(seq_along(ref.tf.align.conv), 10000)

prog.bar.tf <- progress::progress_bar$new(total = length(ref.tf.align.conv[inds.to.explore]))
prog.bar.tf$tick(0)
inds.tf.actions.subset <- lapply(seq_along(ref.tf.align.conv[inds.to.explore]), 
	tf.conv = ref.tf.align.conv[inds.to.explore], tf.rev = ref.tf.align.rev[inds.to.explore],
	ref.pairs.align = ref.pairs.align, ref.actions.db = actions.human.ref.db, 
	function(x, tf.conv, tf.rev, ref.pairs.align, ref.actions.db) {
		one.gp.actions.db <- ref.actions.db[union(which(ref.pairs.align == tf.conv[x]),
			which(ref.pairs.align == tf.rev[x])), ]
		one.gp.gene.ab <- strsplit(tf.conv[x], split = ">", fixed = TRUE)[[1]]
		one.gp.mode.sels <- unique(one.gp.actions.db$mode)
		one.gp.act.sels <- unique(one.gp.actions.db$action)
		ret.orig.inds <- integer()
		for (i in one.gp.mode.sels) {
			for (j in one.gp.act.sels) {
				for (k in 1:2) {  # k = 1, check conv, k = 2, check rev
					inds.chk.exist <- Reduce(intersect ,list(which(one.gp.actions.db$mode == i), 
						which(one.gp.actions.db$action == j),
						which(one.gp.actions.db$is_directional == 't'),
						which(one.gp.actions.db$a_is_acting == 'f')))
					if (length(inds.chk.exist) == 0) {
						next
					}
					inside.gp.gene.ab <- one.gp.gene.ab
					if (k == 2) {
						inside.gp.gene.ab <- c(inside.gp.gene.ab[2], inside.gp.gene.ab[1])
					}
					inds.it <- intersect(which(one.gp.actions.db$inter.GeneName.A == inside.gp.gene.ab[1]),
							which(one.gp.actions.db$inter.GeneName.B == inside.gp.gene.ab[2]))
					inds.it <- intersect(inds.it, inds.chk.exist)
					one.gp.sel.sub.db <- one.gp.actions.db[inds.it, ]
					if (length(inds.it) == 0) {  # rev may have no t-f situation
						# doing nothing
					} else {
						inds.topsel <- intersect(which(one.gp.sel.sub.db$is_directional == 't'), which(one.gp.sel.sub.db$a_is_acting == 't'))
						inds.secsel <- intersect(which(one.gp.sel.sub.db$is_directional == 't'), which(one.gp.sel.sub.db$a_is_acting == 'f'))
						inds.lastsel <- intersect(which(one.gp.sel.sub.db$is_directional == 'f'), which(one.gp.sel.sub.db$a_is_acting == 'f'))
						inds.error.chk <- intersect(which(one.gp.sel.sub.db$is_directional == 'f'), which(one.gp.sel.sub.db$a_is_acting == 't'))
						if (length(inds.error.chk) > 0 || length(inds.topsel) > 0 || length(inds.lastsel) > 0) {
							stop("Unexpected database items after select t-f!")
						}
						inds.sel.ret <- integer(0)
						if (length(inds.secsel) > 0) {
							# check if the rev has the topsel for the same mode and action
							tmp.inds.it.rev <- Reduce(intersect, list(which(one.gp.actions.db$inter.GeneName.A == inside.gp.gene.ab[2]), 
								which(one.gp.actions.db$inter.GeneName.B == inside.gp.gene.ab[1]),
								which(one.gp.actions.db$is_directional == 't'),
								which(one.gp.actions.db$a_is_acting == 't')))
							tmp.inds.it.rev <- Reduce(intersect, list(tmp.inds.it.rev,
								which(one.gp.actions.db$mode == i), 
								which(one.gp.actions.db$action == j)))  # add that mode and action the same
							if (length(tmp.inds.it.rev) == 0) {  # rev has no explicit definitions, so keep this records
								inds.sel.ret <- inds.secsel
							}
						}
						if (length(inds.sel.ret) > 1) {
							message(paste0(tf.conv[x], " in mode: ", i, ", action: ", j, ", has multiple score!"))
							tmp.db <- one.gp.sel.sub.db[inds.sel.ret, ]
							tmp.db <- tmp.db[order(tmp.db$score, decreasing = TRUE), ]
							ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(tmp.db[1])))
						} else {
							if (length(inds.sel.ret) > 0) {
								ret.orig.inds <- c(ret.orig.inds, as.integer(rownames(one.gp.sel.sub.db[inds.sel.ret, ])))
							}
						}
					}
				}
			}
		}
		prog.bar.tf$tick()
		ret.orig.inds
	})

