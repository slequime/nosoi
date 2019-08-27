#This file holds all the function directly related to the nosoiSim function, including the object constructors

#' Top-level function to use Nosoi
#'
#' @param type specifies which type of pathogen we are interested in, either "single" or "dual"-host (e.g. arboviruses).
#' @param popStructure specifies if the population in which the transmission is to occur is structured ("none", "discrete" or "continuous").
#' @param ... arguments to be passed on to the simulator (see below).
#'
#' @details This function determines which general settings the user wants to use for his simulation.
#' @details All other arguments are passed down to the chosen simulator itself, such as \code{\link{singleNone}}, \code{\link{singleDiscrete}}, \code{\link{singleContinuous}}, \code{\link{dualNone}}, \code{\link{dualDiscrete}} or \code{\link{dualContinuous}}.
#' @export nosoiSim
#' @import data.table
#' @import methods
#' @import stats

nosoiSim <- function(type="single", popStructure="none", ...){

  #Sanity checks -------------
  if (! type %in% c("single","dual")) stop("Type of transmission should be 'single' or 'dual'-host.")
  if (! popStructure %in% c("none","discrete","continuous")) stop("Unrecognized parameters for population structure, should be 'none,'discrete' or 'continuous'.")

  #Loading correct script ------------------
  if (type=="single" && popStructure=="none") {
    output <- singleNone(...)
  }
  if (type=="single" && popStructure=="discrete") {
    output <- singleDiscrete(...)
  }
  if (type=="single" && popStructure=="continuous") {
    output <- singleContinuous(...)
  }

  if (type=="dual" && popStructure=="none") {
    output <- dualNone(...)
  }
  if (type=="dual" && popStructure=="discrete") {
    output <- dualDiscrete(...)
  }

  if(type=="dual" && popStructure=="continuous") {
    output <- dualContinuous(...)
  }

  return(output)
}

#' @title nosoiSim Constructor
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param pres.time current time of the simulation
#' @param type population structure (one of "single or "dual)
#' @param pop.A an object of class \code{nosoiSimOne} for population A
#' @param pop.B an object of class \code{nosoiSimOne} for population B
#'
#'
#' @return An object of class \code{nosoiSim}
#'
#' @keywords internal
##

nosoiSimConstructor <- function(total.time,
                                type = c("single", "dual"),
                                pop.A,
                                pop.B = NULL) {

  type <- match.arg(type)

  res <- list(total.time = total.time,
              type = type,
              host.info.A = pop.A,
              host.info.B = pop.B)

  class(res) <- "nosoiSim"

  return(res)

}

#' @title nosoiSimOne Constructor
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param N.infected number of infected hosts
#' @param table.hosts data.table of hosts
#' @param table.state data.table of hosts movement
#' @param popStructure geographical structure (one of "none, "discrete" or "continuous")
#'
#' @return An object of class \code{nosoiSimOne}
#'
#' @keywords internal
##
nosoiSimOneConstructor <- function(N.infected, table.hosts, table.state, prefix.host,
                                   popStructure = c("none", "discrete", "continuous")) {

  popStructure <- match.arg(popStructure)

  res <- list(N.infected = N.infected,
              table.hosts = table.hosts,
              table.state = table.state,
              prefix.host = prefix.host,
              popStructure = popStructure)

  class(res) <- "nosoiSimOne"

  return(res)

}
