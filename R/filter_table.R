#' Table rows filtering
#'
#' Return a logical vector for filtering the table's rows following these steps:
#' \enumerate{
#'   \item Calculate summary values for the rows, either for the whole row or
#'   for each group;
#'   \item Apply 'test' to the summarized values, comparing it to reference
#'   value(s) or using a custom function;
#'   \item If \code{test} was applied by group, reduce the logical result by
#'   row-wise \code{\link[base]{all}}, \code{\link[base]{any}} or
#'   \code{\link[matrixStats]{count}}.
#' }
#'
#' This function can be called in two different ways.  The first and simpler one
#' is using an expression:
#'
#' \code{filter_table(x, mean > 0) # TRUE if the row mean is greater than 0}
#'
#' The second, more flexible form, is using the individual parameters:
#'
#' \code{filter_table(x, summary='mean', test='>', 0) # argument names can be ommited}
#'
#' @param summary  An expression with the summary method and the logical
#' comparison, or the name of a summary method.  See \strong{Details}.  Defaults
#' to \code{"sum"}.
#' @param test  A single character string representing a logical comparison, or
#' a funtion that will be called with the summary vector as the first argument.
#' Defaults to \code{"=="}.
#' @param ...  A value to be compared to the summarized variables, or further
#' potential arguments passed to \code{test} function.  If testing by group, it
#' should be a single value or a vector of the same length as the number of
#' groups.  In this case also, \code{test} is applied using
#' \code{\link[base]{mapply}}, so extra arguments to a custom function may be
#' passed through \code{MoreArgs}.
#' @param MoreArgs  A list of arguments passed to \code{test} by
#' \code{\link[base]{mapply}}.
#' @param method  Method to obtain logical test results by row.  Defaults to
#' \code{"row"}.
#' @param group  A ‘factor’ in the sense that \code{as.factor(group)} defines
#' the grouping, or a list of such factors in which case their
#' \code{\link[base]{interaction}} is used for the grouping.
#' @param drop  Character vector of column names to be ignored.
#' @inheritParams format_table
#' @inheritParams matrixStats::rowMedians
#'
#' @return  A logical vector indicating if each row passed the test.
#'
#' @examples
#' data(geneData, geneCovariate, package='Biobase')
#' gData <- head(geneData, n=100)
#'
#' # Select if standard deviation of log expression level is greater than 0.5.
#' variant <- filter_table(log2(gData), sd > 0.5, na.rm=TRUE)
#'
#' # Consider expressed if the median value is positive in at least one group.
#' expressed <- filter_table(gData, median > 0, method='any',
#'                           group=geneCovariate[c('sex', 'type')])
#'
#' filtered <- gData[variant & expressed, ]
#' nrow(filtered)
#'
#' @export

filter_table <- function(x, summary=c('mean', 'median', 'sd', 'max', 'min',
                                      'sum', 'prod', 'all', 'any', 'count'),
                         test=c('==', '!=', '>', '<', '>=', '<='), ..., MoreArgs=NULL,
                         method=c('row', 'all', 'any', 'count'), group=NULL,
                         na.rm=FALSE, drop=NULL) {

    # Calling form using an expression.
    if (missing(test)) {
        expr <- deparse(substitute(summary))  # summary, comparison operator and value
        args <- regmatches(expr, regexec('(\\w+) ([=!><]=?) (.+)', expr))[[1L]]
        if (!length(args)) stop("Invalid expresssion: ", expr)
        if (length(list(...))) warning("Some arguments were ignored!")
        comp_val <- eval.parent(parse(text=args[4L]))
        return(filter_table(x, args[2L], args[3L], comp_val, 
                            method=method, group=group, na.rm=na.rm, drop=drop))
    }

    # Get necessary functions.
    summary <- switch(match.arg(summary), mean=rowMeans, sum=rowSums,
                      median=matrixStats::rowMedians, sd=matrixStats::rowSds,
                      max=matrixStats::rowMaxs, min=matrixStats::rowMins,
                      prod=matrixStats::rowProds, all=matrixStats::rowAlls,
                      any=matrixStats::rowAnys, count=matrixStats::rowCounts)
    summarize <- function(x) summary(as.matrix(x), na.rm=na.rm)
    if (!is.function(test)) test <- match.fun(match.arg(test))
    method <- switch(match.arg(method), row=NULL, all=matrixStats::rowAlls,
                     any=matrixStats::rowAnys, count=matrixStats::rowCounts)

    # Remove key columns and any others specified by 'drop'.
    if (inherits(x, 'data.table')) drop <- union(drop, data.table::key(x))
    drop <- match(drop, colnames(x), nomatch=0)
    if (any(drop)) x <- `[.data.frame`(x, -drop)

    # Summarize variables and apply test.
    if (is.null(method)) {
        test(summarize(x), ...)
    } else {
        summ <- split_apply(x, FUN=summarize, FACTOR=group)
        summ <- as.matrix(mapply(test, summ, ..., MoreArgs=MoreArgs))
        # If 'x' is a 1-row table, as.matrix(mapply(x)) returns a 1-col matrix.
        if (nrow(x) == 1) summ <- t(summ)
        method(summ)
    }
}
