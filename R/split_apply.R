#' Apply a function to groups of columns
#'
#' Divide the columns of a data frame by the groups defined in \code{FACTOR},
#' apply the function \code{FUN} to each subgroup and then join the resulting
#' vectors in a new data frame.
#'
#' @param X  An object coercible to data frame.
#' @param FUN  A function that receives a data frame and returns a vector.
#' @param FACTOR  A ‘factor’ in the sense that \code{as.factor(FACTOR)} defines
#' the grouping, or a list of such factors in which case their
#' \code{\link[base]{interaction}} is used for the grouping.
#' @param SEP  String to construct the new level labels by joining the
#' constituent ones.  Defaults to \code{"."}.
#' @param ...  Optional arguments to \code{FUN}.
#'
#' @return A data frame, result of the split-apply-combine process.
#'
#' @seealso \code{\link[base]{split}}, \code{\link[base]{lapply}}
#'
#' @examples
#' data(geneData, geneCovariate, package='Biobase')
#' gData <- head(geneData, n=10)
#'
#' meanExprs <- split_apply(log2(gData), FUN=rowMeans,
#'                          FACTOR=geneCovariate[c('sex','type')], na.rm=TRUE)
#' summary(meanExprs)
#'
#' @export

split_apply <- function(X, FUN, FACTOR, SEP='.', ...) {
    x <- split(unclass(as.data.frame(X)), FACTOR, SEP, drop=TRUE)
    cn <- names(x)  # preserve 'SEP' in names
    x <- lapply(lapply(x, as.data.frame), FUN, ...)
    x <- as.data.frame(x, row.names=row.names(X))
    colnames(x) <- cn
    x
}
