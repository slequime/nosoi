#' @title Get Exiting or Moving individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSimOne}.
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
#' @param res an object of class \code{nosoiSimOne}.
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

#' @title Update table state with exiting individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @param exiting a vector (TRUE/FALSE) of exiting individuals
#' @param pres.time current time
#'
#' @return \code{nosoiSim} object
#'
#' @keywords internal
##

updateTableState <- function(res, exiting, pres.time) {

  exiting.ID <- res$table.hosts[exiting]$hosts.ID
  exiting.state <- res$table.state[,is.na(time.to) & hosts.ID %in% exiting.ID]

  res$table.state[exiting.state, `:=` (time.to = as.numeric(pres.time))]

  return(res)
}
