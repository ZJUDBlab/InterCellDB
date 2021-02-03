

#' Tool for showing gragh
#' 
#' @description 
#' This function is a simple interface to show the graph in results that are 
#' from functions \code{GetResult.*}.
#'
#' @param formatted.result List. It must have one item named "plot".
#' @param ... other plot functions that directly passed using \code{+} in definition of \pkg{ggplot2}.
#' It will make slight changes on the final plotting result.
#'
#' 
#'
#' @import grid grid.newpage grid.draw
#'
#' @export
#'
Tool.ShowGraph <- function(
	formatted.result,
	...
) {
	if (!is.null(formatted.result$plot)) {
		plot.res <- formatted.result$plot
		for (each.item in list(...)) {
			plot.res <- plot.res + each.item
		}
		plot.res  # return
	} else {
		plot.res <- formatted.result$grid.plot
		grid.newpage()
		grid.draw(plot.res)
	}
}





#' Tool for write tables in csv
#' 
#' @description 
#' This function is a simple interface to write the tables in results that are 
#' from functions \code{GetResult.*}
#'
#' @param formatted.result List. It must have one item named "table".
#' @param dir.path Character. The written directory path in operating system.
#' @param ... Other params that will be directly passed to write.csv, use \code{?write.csv} for 
#' further help.
#'
#'
#' @importFrom utils write.csv
#'
#' @export
#'
Tool.WriteTables <- function(
	formatted.result,
	dir.path = ".",
	...
) {
	wtables <- formatted.result$table
	if (length(wtables) == 0) {
		return("No table needs to be written!")
	}
	# format dir.path
	dir.path.end <- substring(dir.path, length(dir.path))
	if (!(dir.path.end == "/" || dir.path.end == "\\")) {  # '/' for mac, etc. '\' for windows, etc.
		dir.path <- paste0(dir.path, "/")  # work effective for most operating system
	}
	for (i in 1:length(wtables)) {
		this.table <- wtables[[i]]
		this.tname <- names(wtables)[i]
		this.tnamelist <- strsplit(this.tname, split = ".", fixed = TRUE)[[1]]
		this.tname.final <- paste0(this.tnamelist, collapse = "_")
		write.csv(this.table, file = paste0(dir.path, this.tname.final, ".csv"), ...)
	}
}





#' split character to be data.frame
#' 
#' @description
#' This function generates n-column data frame by spliting character by some specific letters or phrases.
#'
#' @param to.splits.string Character. The string to be splited. 
#' @param to.split.by Character(1). The string will be splited by this parameter, and it will be directly 
#' passed to parameter \code{split} of function \code{strsplit}. 
#' @param res.colnames Character(2). It gives 2 column names in the result table.
#' The splited string gets to have 2 parts of equal length, and 
#' they will be reconstructed to form a table so that 2 column names are needed. 
#'
#'
#'
#' @export
#'
Tool.SplitToGenDataFrame <- function(
  to.splits.string, 
  to.split.by, 
  res.colnames
) {
	tmp.splits <- strsplit(to.splits.string, split = to.split.by, fixed = TRUE)
	tmp.len.splits <- as.integer(unlist(lapply(tmp.splits, FUN = length)))
	tmp.len.splits <- unique(tmp.len.splits)
	if (length(tmp.len.splits) != 1) {
		stop("Given data cannot be uniformly splited, with different splits: ", 
			paste0(tmp.len.splits, collapse = ", "), ".")
	}
	tmp.merge.elems <- as.character(unlist(tmp.splits))
	tmp.merge.base.index <- seq_len(length(tmp.merge.elems) / tmp.len.splits)
	res.prep.list <- list()
	for (i in seq_len(tmp.len.splits)) {
		tmp.indices <- tmp.merge.base.index * tmp.len.splits - (tmp.len.splits - i)
		res.prep.list <- c(res.prep.list, list(tmp.merge.elems[tmp.indices]))
	}
	res.df <- data.frame(res.prep.list, stringsAsFactors = FALSE)
	if (length(res.colnames) == ncol(res.df)) {
		colnames(res.df) <- res.colnames
	} else {
		warning("Given colnames are not matched with the result columns, and unexpected errors may happen!")
		if (length(res.colnames) > ncol(res.df)) {
			colnames(res.df) <- res.colnames[seq_len(ncol(res.df))]
		} else {
			if (length(res.colnames) > 0) {
				colnames(res.df)[seq_along(res.colnames)] <- res.colnames
			}
		}
	}
	return(res.df)
}





# [inside usage]
# get from ?toupper, .simpleCap
Tc.Cap.simple <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
}





# [inside usage]
# parallel version of Tc.Cap.simple
Tc.Cap.simple.vec <- function(to.cap.vec) {
	unlist(lapply(to.cap.vec, FUN = Tc.Cap.simple))
}





#' Speed-up doing unique for data.frame
#'
#' @description
#' This function uses the properties of \code{data.frame} and \code{rownames}, and 
#' makes partial unique process fast and effective.
#'
#' @param xxpairs Data.frame. Any data.frame object.
#' @param cols.select Integer. Specify some columns in integer vector, e.g. c(1:2).
#' Its length \bold{must be >= 2}, or the function will not work properly. 
#'
#' @details
#' When encountering multi-columns tables, e.g. 6 columns or more, \code{unique()} will 
#' be really slow if it is applied on the whole table. However, in most circumstances, 
#' there is no need to apply \code{unique()} on all columns, i.e. only do unique on some 
#' columns, which is exactly the thing this function does.
#'
#'
#'
#' @export
#'
DoPartUnique <- function(
	xxpairs,
	cols.select = c(1:2)
) {
	if (sum(cols.select %in% c(1:ncol(xxpairs))) != length(cols.select)) {
		stop("Columns selected are undefined! Please check again!")
	}
	# rownames(xxpairs) <- NULL
	tmp.uni <- xxpairs[, cols.select]
	rownames(tmp.uni) <- NULL
	tmp.uni <- unique(tmp.uni)
	xxpairs <- xxpairs[as.integer(rownames(tmp.uni)),]
	xxpairs
}



#' Permutate reverse order of odds and evens
#'
#' @description
#' This function generates reverse premutation of odd values and even values 
#' in continues values, e.g. \code{1:10}.
#'
#' @param ncols Integer. An arbitrary integer.
#'
#' @examples
#' # for 1:10 
#' ReverseOddEvenCols(10)
#' # the reuslt is c(2,1,4,3,6,5,8,7,10,9)
#'
#'
#'
#' @export
#'
ReverseOddEvenCols <- function(
	ncols
) {
	len.all <- ncols %/% 2
	val.even <- 2 * 1:len.all
	val.odds <- (2 * 1:len.all) -1
	new.serial <- NULL
	for (i in 1:len.all) {
		new.serial <- c(new.serial, val.even[i], val.odds[i])
	}
	if (ncols %% 2 != 0) {  # odds cols
		# leave last one unchanged
		new.serial <- c(new.serial, ncols)
	}
	#end# return
	new.serial
}





