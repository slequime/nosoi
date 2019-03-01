#' @title Number of Infected units at time t
#'
#' @description
#' For a given time t, this function returns the number of infected active units.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#' @param t time (integer)
#'
#' @return Number of infected units at time t
#'
#' @seealso \code{\link{nosoiSim}}
##
numberInfected <- function(table.nosoi, t) {
  # Attention: need to test that t < t_max
  already_infected <- table.nosoi$inf.time < t  # Strict  inequality: if infected at time t, becomes active at time t + 1
  still_infected <- table.nosoi$out.time > t    # Strict inequality: if out at time t, not active for generation t
  still_infected[is.na(still_infected)] <- TRUE
  return(sum(already_infected & still_infected))
}

#' @title Number of Infected units at time t (BGW)
#'
#' @description
#' For a given time t, this function returns the number of infected active units.
#' The difference with \code{\link{numberInfected}} is that it counts an "out"
#' individual as still there, but with no children.
#' This is for comparision with the BGW process, should be internal only.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#' @param t time (integer)
#'
#' @return Number of infected units at time t
#'
#' @seealso \code{\link{nosoiSim}}
#'
#' @keywords internal
##
numberInfectedBGW <- function(table.nosoi, t) {
  already_infected <- table.nosoi$inf.time <= t  # Non strict inequality: if infected at time t, counts as a member of generation t
  still_infected <- table.nosoi$out.time > t    # Strict inequality: if out at time t, not active for generation t
  still_infected[is.na(still_infected)] <- TRUE
  return(sum(already_infected & still_infected))
}
