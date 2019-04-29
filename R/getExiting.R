#' @title Get Exiting individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @param pExitParsed parsed exit function
#'
#' @return Boolean vector
#'
#' @keywords internal
##
getExiting <- function(res, pres.time, pExitParsed) {

  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

  if (any(active.hosts)) {

    fun <- function(z) {
      pExitParsed$vect(prestime = pres.time, z[, pExitParsed$vectArgs, with = FALSE])
    }
    p.exit.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][["V1"]]

    exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
  }

  exiting.full <- active.hosts
  exiting.full[exiting.full] <- exiting

  return(exiting.full)
}
