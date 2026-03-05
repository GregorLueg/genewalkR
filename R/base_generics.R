# shared generics --------------------------------------------------------------

#' Get the statistical results
#'
#' @param object The `GeneWalk` (see [GeneWalk()]) or `EmbedDrift` class (see
#' [EmbedDrift()]).
#'
#' @returns If found, the statistics results within the class
#'
#' @export
get_stats <- S7::new_generic(
  name = "get_stats",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)
