# This file documents all datasets in this pacakge.
# written in 2020.03.28

### --- Human ---


#' Reference database of human genes
#'
#' @description
#' A dataset contains the basic informations of 41745 human genes.
#' 
#' format List.
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
#' A dataset contains 79373 rows of subcellular locations, which covers 17736 human genes.
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
#' A dataset contains the experiments-based interaction pairs of human genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 214368 rows of pairs in total.
#'
#' @format Data.frame.
#' \itemize{
#'   \item inter.GeneID.*: the unique geneID for every gene in NCBI gene database.
#'    \item inter.GeneName.*: the authorized name for every gene in NCBI gene database.
#'    \item inter.Experiments.Score: score by experiments, which is the <experiments> score in STRING database.
#'    \item inter.Database.Score: score by curated database, which is the <database> score in STRING database.
#'    \item inter.Predicted.Score: score by prediction, which is the columns other than <experiments> and <database.>, 
#'          e.g. gene fusion, gene neighbor, etc.
#' }
#'
#' @details
#' The subset of STRING data used are <taxonomy = Homo Sapiens>, <score[experiments] > 0>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.human.experiments.db"


#' Knowledge database of interaction pairs of human genes
#'
#' @description
#' A dataset contains the curated-database-based interaction pairs of human genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 264689 rows of pairs in total.
#'
#' @format See @format in \code{?pairs.human.experiments.db}.
#'
#' @details
#' The subset of STRING data used are <taxonomy = Homo Sapiens>, <score[combined_score] > 699>, 
#' <score[experiments] = 0>, <score[database] > 699>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.human.knowledge.db"


#' Prediction database of interaction pairs of human genes
#'
#' @description
#' A dataset contains the prediction-based interaction pairs of human genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 73379 rows of pairs in total.
#'
#' @format See @format in \code{?pairs.human.experiments.db}.
#'
#' @details
#' The subset of STRING data used are <taxonomy = Homo Sapiens>, <score[combined_score] > 699>, 
#' <score[experiments] = 0>, <score[database] < 700>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.human.prediction.db"


#' Actions database for human
#'
#' @description
#' A dataset collects detailed action informations about gene-gene pairs. It is extracted from
#' STRING database, and has 3273046 in total.
#'
#' @source STRING <https://string-db.org>
#'
"actions.human.ref.db"

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
#' A dataset contains 66719 rows of subcellular locations, which covers 15605 mouse genes.
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


#' Experiments database of interaction pairs of mouse genes
#'
#' @description
#' A dataset contains the experiments-based interaction pairs of mouse genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 2344 rows of pairs in total.
#'
#' @format Data.frame.
#' \itemize{
#'   \item inter.GeneID.*: the unique geneID for every gene in NCBI gene database.
#'    \item inter.GeneName.*: the authorized name for every gene in NCBI gene database.
#'    \item inter.Experiments.Score: score by experiments, which is the <experiments> score in STRING database.
#'    \item inter.Database.Score: score by curated database, which is the <database> score in STRING database.
#'    \item inter.Predicted.Score: score by prediction, which is the columns other than <experiments> and <database.>, 
#'          e.g. gene fusion, gene neighbor, etc.
#' }
#'
#' @details
#' The subset of STRING data used are <taxonomy = Mus Musculus>, <score[experiments] > 0>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.mouse.experiments.db"


#' Knowledge database of interaction pairs of mouse genes
#'
#' @description
#' A dataset contains the curated-database-based interaction pairs of mouse genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 42146 rows of pairs in total.
#'
#' @format See @format in \code{?pairs.mouse.experiments.db}.
#'
#' @details
#' The subset of STRING data used are <taxonomy = Mus Musculus>, <score[combined_score] > 699>, 
#' <score[experiments] = 0>, <score[database] > 699>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.mouse.knowledge.db"


#' Prediction database of interaction pairs of mouse genes
#'
#' @description
#' A dataset contains the prediction-based interaction pairs of mouse genes, 
#' and is extracted from STRING database. It is further matched by Ligand/Receptor/ECM databases
#' in this package, and as a result, it has 7557 rows of pairs in total.
#'
#' @format See @format in \code{?pairs.mouse.experiments.db}.
#'
#' @details
#' The subset of STRING data used are <taxonomy = Mus Musculus>, <score[combined_score] > 699>, 
#' <score[experiments] = 0>, <score[database] < 700>.
#'
#' @source STRING <https://string-db.org>, HGNC <https://www.genenames.org>, GeneCards <https://www.genecards.org>, etc.
#'
"pairs.mouse.prediction.db"


#' Actions database for mouse
#'
#' @description
#' A dataset collects detailed action informations about gene-gene pairs. It is extracted from
#' STRING database, and has 496206 in total.
#'
#' @source STRING <https://string-db.org>
#'
"actions.mouse.ref.db"


