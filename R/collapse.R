#' Median of a numeric vector
#'
#' Calculates the median without checking the input or removing NAs.
#'
#' @param x
#'
#' @return The median.
#'
#' @seealso \code{\link[stats]{median}}

fmedian <- function(x) {
    # Borrowed from stats::median, without the overhead.
    n <- length(x)
    if (n == 1L) x
    else {
        half <- (n + 1L) %/% 2L
        if (n %% 2L == 1L) sort(x, partial = half)[half]
        else               mean(sort(x, partial = half + 0L:1L)[half + 0L:1L])
    }
}


#' Collapse technical replicates to unique values
#'
#' Summarize replicated data to unique values per identifier using the the
#' median or the mean values.
#'
#' @param id_col  Column name of identifiers column (key column).
#' @param replicates  A vector of length \code{nrow(x)} of replicates
#' identifiers.  Optional, if \code{id_col} isn't in \code{x} or to substitute
#' \code{id_col}.
#' @param method  Method to summarize rows.  Defaults to \code{"median"}.
#' @param drop  Columns to drop (before summarizing and from result).
#' @inheritParams format_table
#'
#' @return Table with unique identifiers and single values for each variable and
#' observation.
#'
#' @note Uses median without checking or data.tables' internaly optimized mean.
#'
#' @seealso \code{\link{collapse_genes}}
#'
#' @import data.table
#' @export

collapse_replicates <- function(x, id_col=NULL, replicates=NULL, method=c('mean', 'median'),
                                drop=character(), df=FALSE) {

    method <- match.arg(method)

    x <- .format_table(x, key_col=id_col, key_data=replicates, default_name='probe')
    ret <- switch(method,
        mean =   x[, lapply(.SD, mean),    by=eval(key(x)[1L]), .SDcols=-drop],
        median = x[, lapply(.SD, fmedian), by=eval(key(x)[1L]), .SDcols=-drop]
    )
    if (df) setDF(ret) else ret
}


#' Select probe sets to represent genes
#'
#' Select a "probe set" to represent each "gene" choosing those with highest or
#' lowest mean expression level (see details in \code{method}).
#'
#' @param genes  A vector of length \code{nrow(x)} of the "genes" identifiers or
#' a mapping of probes to genes (an object coercible to a data table with
#' "probes" identifiers in the first column and "genes" identifiers in the
#' second).  In this last case, \code{x} should have its "probes" identifiers in
#' the \code{id_col} column or in its row names.
#' @param method  Method to choose the "probe set" representative of a "gene":
#' the one with maximum row mean ('maxMean') or with minimum row mean
#' ('minMean'). Defaults to \code{"maxMean"}.
#' @inheritParams format_table
#' @inheritParams collapse_replicates
#'
#' @return Table with unique observations for each "gene".
#'
#' @seealso \code{\link{collapse_replicates}}
#'
#' @import data.table
#' @export

collapse_genes <- function(x, id_col=NULL, genes=NULL,
                           method=c('maxMean', 'minMean'), drop=character(), df=FALSE) {

    selection_method <- switch(match.arg(method), maxMean='last', minMean='first')

    # Check if 'genes' is a mapping: probe -> gene.
    if (ncol(map <- as.data.table(genes)) == 2) {
        x <- .format_table(x, key_col=id_col, default_name='probe')
        setkeyv(map, names(map)[1L])
        genes <- map[x[, 1L, with=F]][[2L]]
    }

    x <- .format_table(x, key_col=id_col, key_data=genes, default_name='gene', sort=FALSE)
    id_col <- colnames(x)[1L]

    # Select probes representative of genes after choosed method.
    x[, rowmean := rowMeans(.SD), .SDcols=-c(id_col, drop)]
    setkey(x, rowmean)
    x[, rowmean := NULL]
    setkeyv(x, id_col)
    ret <- x[unique(x[[id_col]]), mult=selection_method]
    if (df) setDF(ret) else ret
}
