#' Transform dataset to a standard table format
#'
#' Transforms a matrix, data frame or data table to an ordered data table with
#' the key as the first column.  The key column can be passed, selected by its
#' column name or extracted from row names.
#'
#' @param x  A matrix, data frame or data table.
#' @param key_col  Character vector of one or more column names which is passed
#' to \code{\link[data.table]{setkey}}.  If \code{x} has a key, it can be
#' ommited.
#' @param key_data  A vector of length \code{nrow(x)} to be used as key.  If the
#' (first) column name in \code{key_col} is already a column of \code{x}, it
#' will be substituted by \code{key_data} and the table reordered.  Otherwise,
#' and if \code{key_data} is  missing, \code{row.names(x)} will be used instead.
#' @param sort  Wheter to sort by \code{key_col}.  Defaults to \code{TRUE}.
#' @param df  Return a data frame instead of a data table.  Defaults to
#' \code{FALSE}.
#'
#' @return A data table (or data frame) based on \code{x} and ordered by \code{key_col}.
#'
#' @examples
#' data(geneData, geneCovariate, package='Biobase')
#' gData <- head(geneData, n=5)
#'
#' gData
#' class(gData)
#'
#' gData <- format_table(gData, 'Probeset')
#' gData
#' class(gData)
#' key(gData)
#'
#' gData <- format_table(gData, key_data=row.names(geneData)[6:10], sort=FALSE)
#' gData
#' key(gData)
#'
#' @import data.table
#' @export

format_table <- function(x, key_col=key(x), key_data=NULL, sort=TRUE, df=FALSE)
{ .format_table(x, key_col, key_data, default_name=NULL, sort, df) }

.format_table <- function(x, key_col=key(x), key_data=NULL, default_name=NULL, sort=TRUE, df=FALSE) {

    # Guarantee a 'key_col'.
    if (is.null(key_col)) {
        if (is.null(default_name)) stop("Must have a 'key_col'!")
        key_col <- default_name
    }
    # Transform 'x' to data.table.
    if (!is.data.table(x)) {
        # Do we need to extract the key column from row.names?
        keep_rn <- is.null(key_data) && !key_col[1L] %in% colnames(x)
        if (keep_rn) {
            # Does 'x' have row.names?
            rn <- rownames(x)
            if (is.null(rn) || identical(rn, as.character(seq_along(rn)))) {
                stop("'x' doesn't have row.names to be used as the key column.")
            }
        }
        keep_rn <- if (keep_rn) key_col else FALSE
        x <- switch(class(x),
            data.frame = setDT(x, keep.rownames=keep_rn),
            matrix = as.data.table(x, keep.rownames=keep_rn)
        )
    }

    # Add key column to the table.
    if (!is.null(key_data)) set(x, j=key_col[1L], value=key_data)

    # Sort table by key column(s).
    if (sort && !identical(key(x), key_col)) setkeyv(x, key_col)

    # Set key column(s) in the first position(s).
    key_pos <- match(key_col, names(x))
    setcolorder(x, c(key_col, names(x)[-key_pos]))
    if (df) setDF(x) else x
}
