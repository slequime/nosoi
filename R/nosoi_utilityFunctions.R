#' @title Progress bar
#'
#' @description
#' Echos the state of the simulation at any given time step provided by the user.
#'
#' @param Host.count number of infected hosts
#' @param pres.time current time of the simulation
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
#' @param length.sim the length (in unit of time) over which the simulation should be run.
#' @param max.infected the maximum number of hosts that can be infected in the simulation.
#'
#'
#' @keywords internal
##

progressMessage <- function(Host.count, pres.time, print.step, length.sim, max.infected) {
  if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count," (",round((Host.count/max.infected)*100,digits=0),"% of maximum infected hosts).")}
}

#' @title End message
#'
#' @description
#' Message that ends the simulation
#'
#' @param Host.count number of infected hosts
#' @param pres.time current time of the simulation
#'
#'
#' @keywords internal
##

endMessage <- function(Host.count, pres.time) {
  message("done. \nThe simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")
}

#' @title nosoiSim Constructor
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param N.infected number of infected hosts
#' @param pres.time current time of the simulation
#' @param table.hosts data.table of hosts
#' @param table.state data.table of hosts movement
#' @param type name of the diffusion
#'
#' @return An object of class \code{nosoiSim}
#'
#' @keywords internal
##

nosoiSimConstructor <- function(N.infected, pres.time, table.hosts, table.state,
                                type = c("singleNone", "singleDiscrete", "singleContinuous")) {

  type <- match.arg(type)

  res <- list(total.time = pres.time,
              N.infected = N.infected,
              table.hosts = table.hosts,
              table.state = table.state,
              type = type)

  class(res) <- "nosoiSim"

  return(res)

}

#' @title get Position Infected
#'
#' @description
#' Return the relevent position of the infected individual.
#'
#' @param nosoiSim an object of class \code{nosoiSim}.
#' @param df.meetTransmit current data.table for transmission
#' @param i individual
#'
#' @return A vector of positions
#'
#' @keywords internal
##
getPositionInfected <- function(nosoiSim, df.meetTransmit, i) {
  if (nosoiSim$type == "singleNone") return(NA)
  if (nosoiSim$type == "singleDiscrete") return(df.meetTransmit[i, ]$current.in)
  if (nosoiSim$type == "singleContinuous") return(c(df.meetTransmit[i, ]$current.in.x, df.meetTransmit[i, ]$current.in.y))
}

#' @title Should we build the table.host table
#'
#' @param nosoiSim an object of class \code{nosoiSim}.
#'
#' @return boolean
#'
#' @keywords internal
##
keepState <- function(nosoiSim) {
  if (nosoiSim$type == "singleNone") return(FALSE)
  if (nosoiSim$type == "singleDiscrete") return(TRUE)
  if (nosoiSim$type == "singleContinuous") return(TRUE)
}
