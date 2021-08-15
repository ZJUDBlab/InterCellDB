# This is the R scripts designed for generating more accurate-matching-to-STRING 
# ensembl database.
# 
# The reason is ensembl Protein.stable.ID get from STRING is hard to perfectly match
# to ensembl release 101 - GRCh38 & GRCh37.
# There is some unmatched ones that is considered to be manually checked and got them 
# matched anyway.


###### human
# manual patches, which is superior to release 101 database, as it is checked in 
# ensembl website maually.
# This is the merged result
manual.patch <- readRDS("../ensembl-manual-patches/human-GRch-STRING-patches-upd.rds")
#

# GRCh38
ensembl.r101.g38.9606.db <- read.csv("../raw-db-set/GRCh38p13-release101-human.txt", header = TRUE)
colnames(ensembl.r101.g38.9606.db)[4] <- "NCBI.gene.ID"
ensembl.r101.g38.9606.db <- ensembl.r101.g38.9606.db[which(ensembl.r101.g38.9606.db$Protein.stable.ID != ""), ]
# 2020-12-01, 113141 rows

# GRCh37
ensembl.g37.9606.db <- read.csv("../raw-db-set/GRCh37-human.txt", header = TRUE)
colnames(ensembl.g37.9606.db)[4] <- "NCBI.gene.ID"
ensembl.g37.9606.db <- ensembl.g37.9606.db[which(ensembl.g37.9606.db$Protein.stable.ID != ""), ]
# 2020-12-01, 113750 rows

# GRCh38 removes those Protein.stable.ID in manual patch
m.patch.Protein.id <- levels(factor(manual.patch$old))
inds.t38.in <- which(ensembl.r101.g38.9606.db$Protein.stable.ID %in% m.patch.Protein.id)
ensembl.r101.g38.9606.db <- ensembl.r101.g38.9606.db[setdiff(1:nrow(ensembl.r101.g38.9606.db), inds.t38.in), ]
# 2020-12-01, 113141 rows, not changed here

# GRCh37 removes those Protein.stable.ID in manual patch & GRCh38 r101
inds.t37.in <- which(ensembl.g37.9606.db$Protein.stable.ID %in% m.patch.Protein.id)
ensembl.g37.9606.db <- ensembl.g37.9606.db[setdiff(1:nrow(ensembl.g37.9606.db), inds.t38.in), ]
# 2020-12-01, 113750 rows
r101.g38.Protein.id <- levels(factor(ensembl.r101.g38.9606.db$Protein.stable.ID))
inds.g37.x.g38 <- which(ensembl.g37.9606.db$Protein.stable.ID %in% r101.g38.Protein.id)
ensembl.g37.9606.db <- ensembl.g37.9606.db[setdiff(1:nrow(ensembl.g37.9606.db), inds.g37.x.g38), ]

# merge all 
manual.patch <- manual.patch[, c(2,1,3,4)]
colnames(manual.patch)[1:3] <- colnames(ensembl.r101.g38.9606.db)[1:3]
manual.patch <- manual.patch[which(!is.na(manual.patch$Gene.stable.ID)), ]
manual.patch <- manual.patch[setdiff(1:nrow(manual.patch), which(manual.patch$If.in.entrez == FALSE)), ]
# merge result
ensembl.9606.db <- rbind(ensembl.r101.g38.9606.db[, 1:3], rbind(ensembl.g37.9606.db[, 1:3], manual.patch[, 1:3]))
for (i in 1:3) {
	ensembl.9606.db[, i] <- as.character(ensembl.9606.db[, i])
}
# 2020-12-01, 131631 rows





###### mouse

# get from STRING db, all the ENSG**
mouse.all.string.v11.proteins <- readRDS("../raw-db-set/unmatches-ensembl-string-v11-mouse.rds")

# formatted colnames
std.colnames.GRCm <- c("Gene.stable.ID", "Protein.stable.ID", "Gene.name", "NCBI.gene.ID", "MGI.ID")

# GRCmp6 realease 101 as the base
ensembl.p6.db <- read.csv("../raw-db-set/GRCm38p6-release101-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(ensembl.p6.db) <- std.colnames.GRCm
ensembl.p6.db <- ensembl.p6.db[which(ensembl.p6.db$Protein.stable.ID != ""), ]
# 115173 rows -> 69022 rows

# GRCmp5 release 91 as the backup 1
ensembl.p5.db <- read.csv("../raw-db-set/GRCm38p5-release91-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(ensembl.p5.db) <- std.colnames.GRCm
ensembl.p5.db <- ensembl.p5.db[which(ensembl.p5.db$Protein.stable.ID != ""), ]
# 108738 rows -> 65693 rows

# GRCmp4 release 86 as the backup 2
ensembl.p4.db <- read.csv("../raw-db-set/GRCm38p4-release86-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(ensembl.p4.db) <- std.colnames.GRCm
ensembl.p4.db <- ensembl.p4.db[which(ensembl.p4.db$Protein.stable.ID != ""), ]
# 97912 rows -> 61122 rows

# GRCmp3 release 80 as the backup 3
ensembl.p3.db <- read.csv("../raw-db-set/GRCm38p3-release80-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(ensembl.p3.db) <- std.colnames.GRCm
ensembl.p3.db <- ensembl.p3.db[which(ensembl.p3.db$Protein.stable.ID != ""), ]
# 88998 rows -> 56550 rows

# GRCmp2 release 77 as the backup 4
ensembl.p2.db <- read.csv("../raw-db-set/GRCm38p2-release77-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(ensembl.p2.db) <- std.colnames.GRCm
ensembl.p2.db <- ensembl.p2.db[which(ensembl.p2.db$Protein.stable.ID != ""), ]
# 90368 rows -> 62269 rows

### match process
# 1. p6 r101 matching
inds.match.p6.101 <- which(mouse.all.string.v11.proteins %in% ensembl.p6.db$Protein.stable.ID)
# 2. p5 r91 matching
inds.match.p5.r91 <- which(mouse.all.string.v11.proteins %in% ensembl.p5.db$Protein.stable.ID)
inds.match.p5.r91 <- setdiff(inds.match.p5.r91, inds.match.p6.101)  # p6 get the highest priority
# 3. p4 r91 matching
inds.match.p4.r86 <- which(mouse.all.string.v11.proteins %in% ensembl.p4.db$Protein.stable.ID)
inds.match.p4.r86 <- setdiff(inds.match.p4.r86, inds.match.p6.101)
inds.match.p4.r86 <- setdiff(inds.match.p4.r86, inds.match.p5.r91)
# 4. p3 r80 matching
inds.match.p3.r80 <- which(mouse.all.string.v11.proteins %in% ensembl.p3.db$Protein.stable.ID)
inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p6.101)
inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p5.r91)
inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p4.r86)
# 5. p2 r77 matching
inds.match.p2.r77 <- which(mouse.all.string.v11.proteins %in% ensembl.p2.db$Protein.stable.ID)
inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p6.101)
inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p5.r91)
inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p4.r86)
inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p3.r80)
# merge all matches
inds.match.p2t6 <- c(inds.match.p6.101, inds.match.p5.r91, inds.match.p4.r86, inds.match.p3.r80, inds.match.p2.r77)
#
# [NOTE] as this inds cover all the STRING Protein.stable.ID, so there is no need to do manual patch.
#
mouse.ensembl.patch.list <- list(p5 = ensembl.p5.db[which(ensembl.p5.db$Protein.stable.ID %in% mouse.all.string.v11.proteins[inds.match.p5.r91]), ], 
	p4 = ensembl.p4.db[which(ensembl.p4.db$Protein.stable.ID %in% mouse.all.string.v11.proteins[inds.match.p4.r86]), ],
	p3 = ensembl.p3.db[which(ensembl.p3.db$Protein.stable.ID %in% mouse.all.string.v11.proteins[inds.match.p3.r80]), ],
	p2 = ensembl.p2.db[which(ensembl.p2.db$Protein.stable.ID %in% mouse.all.string.v11.proteins[inds.match.p2.r77]), ])
ensembl.10090.db <- rbind(ensembl.p6.db[which(ensembl.p6.db$Protein.stable.ID %in% mouse.all.string.v11.proteins[inds.match.p6.101]), ], dplyr::bind_rows(mouse.ensembl.patch.list))
# 21837 rows
ensembl.10090.db <- DoPartUnique(ensembl.10090.db, c(1:3))  # 21291 rows


