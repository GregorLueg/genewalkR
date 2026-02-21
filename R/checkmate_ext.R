# checkmate extensions ---------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Check node2vec parameters
#'
#' @description Checkmate extension for checking the node2vec parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkNode2VecParams <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "p",
      "q",
      "walks_per_node",
      "walk_length",
      "num_workers",
      "batch_size",
      "n_epochs",
      "n_negatives",
      "window_size",
      "lr"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  integer_rules <- list(
    "walks_per_node" = "I1",
    "walk_length" = "I1",
    "num_workers" = "I1",
    "batch_size" = "I1",
    "n_epochs" = "I1",
    "n_negatives" = "I1",
    "window_size" = "I1"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(x, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in node2vec parameters is incorrect:",
          "walks_per_node, walk_length, num_workers, batch_size, n_epochs,",
          "n_negatives, and window_size need to be integers."
        ),
        broken_elem
      )
    )
  }

  numeric_rules <- list(
    "p" = "N1",
    "q" = "N1",
    "lr" = "N1"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(x, numeric_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in node2vec parameters is incorrect:",
          "p, q, and lr need to be numeric."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

#' Assert node2vec parameters
#'
#' @description Checkmate extension for asserting the node2vec parameters.
#'
#' @inheritParams checkNode2VecParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertNode2VecParam <- checkmate::makeAssertionFunction(checkNode2VecParams)

## synthetic data parameter ----------------------------------------------------

#' Check GeneWalk data parameters
#'
#' @description Checkmate extension for checking GeneWalkData parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkGeneWalkDataParams <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  integer_params <- c(
    "n_signal_genes",
    "n_noise_genes",
    "n_roots",
    "depth",
    "branching",
    "min_annotations",
    "max_annotations",
    "min_noise_subtrees"
  )
  numeric_params <- c("p_lateral", "p_ppi")
  res <- checkmate::checkNames(
    names(x),
    must.include = c(integer_params, numeric_params)
  )
  if (!isTRUE(res)) {
    return(res)
  }
  for (param in integer_params) {
    res <- checkmate::qtest(x[[param]], "I1")
    if (!isTRUE(res)) {
      return(sprintf("`%s` must be a single integer", param))
    }
  }
  for (param in numeric_params) {
    res <- checkmate::qtest(x[[param]], "N1")
    if (!isTRUE(res)) {
      return(sprintf("`%s` must be a single numeric", param))
    }
  }
  TRUE
}

#' Assert GeneWalk data parameters
#'
#' @description Checkmate extension for asserting GeneWalkData parameters.
#'
#' @inheritParams checkGeneWalkDataParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertGeneWalkDataParams <- checkmate::makeAssertionFunction(
  checkGeneWalkDataParams
)

## interaction check -----------------------------------------------------------

#' Check GeneWalk input data.table
#'
#' @description Checkmate extension for checking GeneWalkDataTable parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
assertGeneWalkDataTable <- checkmate::makeAssertionFunction(
  checkGeneWalkDataTable
)
