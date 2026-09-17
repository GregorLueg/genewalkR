# checkmate extensions ---------------------------------------------------------

## data checks -----------------------------------------------------------------

#' Check GeneWalk input data.table
#'
#' @description Checkmate extension for checking GeneWalkDataTable parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkGeneWalkDataTable <- function(x) {
  res <- checkmate::checkDataTable(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(names(x), must.include = c("from", "to"))
  if (!isTRUE(res)) {
    return(res)
  }
  TRUE
}

#' Assert GeneWalk input data.table
#'
#' @description Checkmate extension for asserting GeneWalkDataTable parameters.
#'
#' @inheritParams checkGeneWalkDataTable
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertGeneWalkDataTable <- checkmate::makeAssertionFunction(
  checkGeneWalkDataTable
)


#' Check GeneWalk graph data.table
#'
#' @description Checkmate extension for checking the graph_dt parameter of a
#' GeneWalk object.
#'
#' @param x The data.table to check.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkGeneWalkGraphDt <- function(x) {
  res <- checkmate::checkDataTable(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(names(x), must.include = c("from", "to", "type"))
  if (!isTRUE(res)) {
    return(res)
  }

  required_types <- c("hierarchy", "interaction", "part_of")
  missing_types <- setdiff(required_types, unique(x$type))
  if (length(missing_types) > 0) {
    return(sprintf(
      "Column 'type' in graph_dt must include: %s. Missing: %s.",
      paste(required_types, collapse = ", "),
      paste(missing_types, collapse = ", ")
    ))
  }

  TRUE
}

#' Assert GeneWalk graph data.table
#'
#' @description Checkmate extension for asserting the graph_dt parameter of a
#' GeneWalk object.
#'
#' @inheritParams checkGeneWalkGraphDt
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertGeneWalkGraphDt <- checkmate::makeAssertionFunction(checkGeneWalkGraphDt)
