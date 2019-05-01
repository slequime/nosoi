#' @title Get Exiting or Moving individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @param pasedFunction parsed exit/moving function
#'
#' @return Boolean vector
#'
#' @keywords internal
##
getExitingMoving <- function(res, pres.time, pasedFunction) {

  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

  if (any(active.hosts)) {

    p.exitMove.values <- applyFunctionToHosts(res, pres.time, pasedFunction, active.hosts)

    exitMove <- drawBernouilli(p.exitMove.values) #Draws K bernouillis with various probability (see function for more detail)
  }

  if(all(active.hosts == FALSE)) exitMove <- FALSE

  exitMove.full <- active.hosts
  exitMove.full[exitMove.full] <- exitMove

  return(exitMove.full)
}

#' @title Apply a function to table.host
#'
#' @description
#' Return a vector of the results of the function
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @param pasedFunction parsed exit/moving function
#' @param active.hosts a boolean vector of active hosts
#'
#' @return result vector
#'
#' @keywords internal
##
applyFunctionToHosts <- function(res, pres.time, pasedFunction, active.hosts) {
  fun <- function(z) {
    pasedFunction$vect(prestime = pres.time, z[, pasedFunction$vectArgs, with = FALSE])
  }
  return(res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][["V1"]])
}
