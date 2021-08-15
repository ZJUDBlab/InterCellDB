# STRING database-full.links explore
#
# !!!
# [deprecated definition at 2020.01.30]
# Database definitions see inside package.
# !!!
# 


# STRING use Ensembl IDs
# ref db: Ensembl database
#
ensembl.10090.db <- read.csv("../STRING/10090-csv-mgi-mart_export.txt", stringsAsFactors = FALSE)
ensembl.10090.db <- ensembl.10090.db[which(ensembl.10090.db$Protein.stable.ID != ""), ]
ensembl.10090.db <- ensembl.10090.db[which(ensembl.10090.db$MGI.ID != ""), ]

# get from other file "10090-interaction-STRING-full-ext-20191109.rds"
interaction.list.string  # high score > 0.7


# Known Interactions
# as defined in STRING, including
# (1) curated databases; (2) experimentally determined
# defined as "Known"
interaction.list.known.string <- interaction.list.string[
	union(which(interaction.list.string$database > 700), 
		  which(interaction.list.string$experiments > 700)
	), ]
# done, with 546272 rows, at 2019.10.26, deprecated
# done, with 546117 rows, at 2019.11.09, 
# the number change is fine, because 1109 does unique() step, while the 1026 don't
# done, with 547446 rows, at 2019.11.15, as mt-DNAs included

# saveRDS, interaction.list.known.string to be sspairs
sspairs

# Predicted Interactions
# as defined in STRING, including
# (1) gene neighborhood; (2) gene fusions; (3) gene co-occurrence
# defined as "Predicted"
interaction.list.predicted.string <- interaction.list.string[
	union(
		union(which(interaction.list.string$neighborhood > 700),
			  which(interaction.list.string$fusion > 700)),
		which(interaction.list.string$cooccurence > 700)
	), ]
# done, with 120 rows, at 2019.10.26


# Other
# as defined in STRING, including
# (1) textmining; (2) co-expression; (3) protein homology
# defined as "Other"
ind.textm.total <- union(which(interaction.list.string$textmining > 700), which(interaction.list.string$textmining_transferred > 700))
ind.coexp.total <- union(which(interaction.list.string$coexpression > 700), which(interaction.list.string$coexpression_transferred > 700))
ind.trans.1 <- union(which(interaction.list.string$neighborhood_transferred > 700), which(interaction.list.string$homology > 700))
ind.trans.2 <- union(which(interaction.list.string$experiments_transferred > 700), which(interaction.list.string$database_transferred > 700))
ind.other.x1 <- union(ind.trans.1, ind.trans.2)
ind.other.x2 <- union(ind.textm.total, ind.coexp.total)
ind.other.total <- union(ind.other.x1, ind.other.x2)
interaction.list.other.string <- interaction.list.string[ind.other.total, ]
# done, with 133430 rows, at 2019.10.26



## Merge "Predicted" and "Other" to be one database
## and renamed as "presumed"
interaction.list.presumed.string <- rbind(interaction.list.predicted.string, interaction.list.other.string)
# remove known interaction pairs within
ind.mixed.known <- union(which(interaction.list.presumed.string$database > 700), which(interaction.list.presumed.string$experiments > 700))
if (length(ind.mixed.known) > 0)
	interaction.list.presumed.string <- interaction.list.presumed.string[-ind.mixed.known, ]


# saveRDS interaction.list.presumed.string, to be aapairs
aapairs
# rows 73198, at 2019.11.09
# rows 73340, at 2019.11.15, as mt-DNAs included





















