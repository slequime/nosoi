## Draw K bernouillis with probas p1, ..., pK

#' @title Draw newly infected
#'
#' @description
#' For each encounter, simulate wether a new individual is infected.
#'
#' @param p vector of size K, giving the probability that each encounter
#' leads to an infection.
#'
#' @return Boolean vector giving the newly infected individuals.
#' @keywords internal
##
drawBernouilli <- function(p) {
  return(runif(length(p), 0, 1) < p)
}
