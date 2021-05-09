
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
#' @importFrom grid grid.newpage grid.draw
#'
#' @export
#'
Tool.ShowGraph <- function(
	formatted.result,
	...
) {
	if (!is.null(formatted.result$plot)) {
		plot.res <- formatted.result$plot
		# for (each.item in list(...)) {
		# 	plot.res <- plot.res + each.item
		# }
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
