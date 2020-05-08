
nLR.db <- read.csv("/Users/jinziyang/Downloads/2015Nature-LR-list.csv",
	stringsAsFactors = FALSE)
nLigand <- unique(nLR.db$Ligand.ApprovedSymbol)
nReceptor <- unique(nLR.db$Receptor.ApprovedSymbol)

nLigand  # 708
nReceptor  # 691

df.ligand <- data.frame(Gene = nLigand, user.lr.type = "Ligand", stringsAsFactors = FALSE)
df.receptor <- data.frame(Gene = nReceptor, user.lr.type = "Receptor", stringsAsFactors = FALSE)

user.lr.db <- rbind(df.ligand, df.receptor)

user.lr.db <- Tool.AddUserRestrictDB(user.lr.db, genes.human.ref.db)

# prepare to be used
user.lr.db
user.lr.colname <- "user.lr.type"

