# This file documents all datasets in this pacakge.
# written in 2020.03.28

globalVariables(c("Uniprot.key.map.list",
	"genes.human.ref.db", "anno.location.human.ref.db", "anno.type.human.ref.db", 
	"go.human.ref.db", "pairs.human.db", "actions.human.ref.db",
	"genes.mouse.ref.db", "anno.location.mouse.ref.db", "anno.type.mouse.ref.db", 
	"go.mouse.ref.db", "pairs.mouse.db", "actions.mouse.ref.db")
)

### Small Tools

#' Merged Uniprot Keywords
#'
#' @description
#' A dataset contains the mapping of original Uniprot Keywords and remapped keys
#' 
#' @format Data.frame.
#' \itemize{
#'    \item {\code{orig.uniprot.keywords}: } {[TODO]}
#'    \item {\code{merged.molecular.function}: } {[TODO]
#'          }
#' }
#'
#' @source <https://www.uniprot.org>, and download the Keywords-[molecular function] of genes.
#'
"Uniprot.key.map.list"


### --- Human ---

#' Reference database of human genes
#'
#' @description
#' A dataset contains the basic informations of 41745 human genes.
#' 
#' @format List.
#' \itemize{
#'    \item {\code{gene.ncbi.db}: } {it is downloaded from NCBI ftp site:"/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz".}
#'    \item {\code{gene.synonyms.db}: } {data.frame, it records the synonyms-genenames pairs, which are generated from \code{gene.ncbi.db}.
#'          It maps gene synonyms to its authorized name and ID.
#'          }
#'    \item {\code{gene.dup.synonyms.db}: } {data.frame, which documents synonyms with multiple matching genes.}
#' }
#'
#' @source <ftp://ftp.ncbi.nih.gov/gene/>, and file path is /DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz.
#'
"genes.human.ref.db"


#' Annotation database for subcellular location of human genes
#'
#' @description
#' A dataset contains 79628 rows of subcellular locations, which covers 17646 human genes.
#'
#' @source <https://compartments.jensenlab.org>, and download the knowledge channel of human databases.
#'
"anno.location.human.ref.db"


#' Annotation database for molecular function of human genes
#'
#' @description
#' A dataset contains 30756 rows of molecular functions, which covers 20195 human genes.
#'
#' @source <https://www.uniprot.org>, and download the Keywords-[molecular function] of genes.
#'
"anno.type.human.ref.db"


#' Experiments database of interaction pairs of human genes
#'
#' @description
#' A dataset contains the all-from-STRING interaction pairs of human genes.
#' In this package, it has 5758309 unique pairs in total.
#'
#' @format Data.frame.
#' \itemize{
#'   \item inter.GeneID.*: the unique geneID for every gene in NCBI gene database.
#'    \item inter.GeneName.*: the authorized name for every gene in NCBI gene database.
#'    \item inter.Experiments.Score: score by experiments, which is the <experiments> score in STRING database.
#'    \item inter.Database.Score: score by curated database, which is the <database> score in STRING database.
#'    \item inter.Predicted.Score: score by prediction, which is the columns other than <experiments> and <database.>, 
#'          e.g. gene fusion, gene neighbor, etc.
#'    \item inter.Combined.Score: score combing all above scores.
#' }
#'
#' @details
#' The subset of STRING data used are <taxonomy = Homo Sapiens>.
#'
#' @source STRING <https://string-db.org>.
#'
"pairs.human.db"


#' Actions database for human
#'
#' @description
#' A dataset collects detailed action informations about gene-gene pairs. It is extracted from
#' STRING database, and has 3409770 in total.
#'
#' @source STRING <https://string-db.org>
#'
"actions.human.ref.db"


#' GO reference database for human
#'
#' @description
#' A dataset collects human genes and their GO annotations. It is extracted from
#' GO annotations of mixed species, and has 322873 rows in total covering 19606 genes.
#'
#' @source <ftp://ftp.ncbi.nih.gov/gene/>, and file path is /DATA/gene2go.gz.
#'
"go.human.ref.db"


### --- Mouse ---

#' Reference database of mouse genes
#'
#' @description
#' A dataset contains the basic informations of 69864 mouse genes.
#'
#' @format List.
#' \itemize{
#'    \item {\code{gene.ncbi.db}: } {it is downloaded from NCBI ftp site:"/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz".}
#'    \item {\code{gene.synonyms.db}: } {data.frame, it records the synonyms-genenames pairs, which are generated from \code{gene.ncbi.db}.
#'          It maps gene synonyms to its authorized name and ID.
#'          }
#'    \item {\code{gene.dup.synonyms.db}: } {data.frame, which documents synonyms with multiple matching genes.}
#' }
#'
#' @source <ftp://ftp.ncbi.nih.gov/gene/>, and file path is /DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz.
#'
"genes.mouse.ref.db"


#' Annotation database for subcellular location of mouse genes
#'
#' @description
#' A dataset contains 66569 rows of subcellular locations, which covers 15566 mouse genes.
#'
#' @source <https://compartments.jensenlab.org>, and download the knowledge channel of mouse databases.
#'
"anno.location.mouse.ref.db"


#' Annotation database for molecular function of mouse genes
#'
#' @description
#' A dataset contains 25831 rows of molecular functions, which covers 16655 mouse genes.
#'
#' @source <https://www.uniprot.org>, and download the Keywords-[molecular function] of genes.
#'
"anno.type.mouse.ref.db"


#' Database of interaction pairs of mouse genes
#'
#' @description
#' A dataset contains the all-from-STRING interaction pairs of mouse genes.
#' In this package, it has 5882115 unique pairs in total.
#'
#' @format Data.frame.
#' \itemize{
#'    \item inter.GeneID.*: the unique geneID for every gene in NCBI gene database.
#'    \item inter.GeneName.*: the authorized name for every gene in NCBI gene database.
#'    \item inter.Experiments.Score: score by experiments, which is the <experiments> score in STRING database.
#'    \item inter.Database.Score: score by curated database, which is the <database> score in STRING database.
#'    \item inter.Predicted.Score: score by prediction, which is the columns other than <experiments> and <database.>, 
#'          e.g. gene fusion, gene neighbor, etc.
#'    \item inter.Combined.Score: score combing all above scores.
#' }
#'
#' @details
#' The subset of STRING data used are <taxonomy = Mus Musculus>.
#'
#' @source STRING <https://string-db.org>.
#'
"pairs.mouse.db"


#' Actions database for mouse
#'
#' @description
#' A dataset collects detailed action informations about gene-gene pairs. It is extracted from
#' STRING database, and has 4763978 in total.
#'
#' @source STRING <https://string-db.org>
#'
"actions.mouse.ref.db"


#' GO reference database for mouse
#'
#' @description
#' A dataset collects mouse genes and their GO annotations. It is extracted from
#' GO annotations of mixed species, and has 367110 rows in total covering 24193 genes.
#'
#' @source <ftp://ftp.ncbi.nih.gov/gene/>, and file path is /DATA/gene2go.gz.
#'
"go.mouse.ref.db"

