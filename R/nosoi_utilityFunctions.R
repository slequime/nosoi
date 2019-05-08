#' @title Progress bar
#'
#' @description
#' Echos the state of the simulation at any given time step provided by the user.
#'
#' @param Host.count.A number of infected hosts of host A.
#' @param Host.count.B number of infected hosts of host B.
#' @param pres.time current time of the simulation
#' @param print.step if print.progress is TRUE, step with which the progress message will be printed.
#' @param length.sim the length (in unit of time) over which the simulation should be run.
#' @param max.infected.A the maximum number of hosts that can be infected in the simulation for host A.
#' @param max.infected.B the maximum number of hosts that can be infected in the simulation for host B.
#' @param type either single/dual host
#'
#' @keywords internal
##

progressMessage <- function(Host.count.A, Host.count.B=NULL, pres.time, print.step, length.sim, max.infected.A, max.infected.B=NULL, type="single") {

  if(type == "single"){
    if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count.A," (",round((Host.count.A/max.infected.A)*100,digits=0),"% of maximum infected hosts).")}
  }

  if(type == "dual"){
    if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: (A) ", Host.count.A," (",round((Host.count.A/max.infected.A)*100,digits=0),"% of maximum infected hosts); (B) ", Host.count.B," (",round((Host.count.B/max.infected.B)*100,digits=0),"% of maximum infected hosts).")}
  }
}

#' @title End message
#'
#' @description
#' Message that ends the simulation
#'
#' @param Host.count.A number of infected hosts (host A)
#' @param Host.count.B number of infected hosts (host B)
#' @param pres.time current time of the simulation
#' @param type either single/dual host
#'
#' @keywords internal
##

endMessage <- function(Host.count.A, Host.count.B=NULL, pres.time, type="single") {
  if(type == "single"){
    message("done. \nThe simulation has run for ",pres.time," units of time and a total of ",Host.count.A," hosts have been infected.")
  }
  if(type == "dual"){
    message("done. \nThe simulation has run for ",pres.time," units of time and a total of ",Host.count.A," (A) and ",Host.count.B, " (B) hosts have been infected.")
  }
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

#' @title get host info
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param res an object of class \code{nosoiSim}
#' @param what the information to get. One of "table.hosts", "N.infected", "table.state", "popStructure"
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with host informations
#'
#' @export
#'
##
getHostInfo <- function(res,
                        what = c("table.hosts", "N.infected", "table.state", "popStructure"),
                        pop = "A") {

  what <- match.arg(what)

  if (pop == "A") return(res$host.info.A[[what]])
  if (pop == "B") return(res$host.info.B[[what]])

  stop(paste0("Population ", pop, " is not recognized"))
}

#' @title Get table hosts
#'
#' @description
#' Get "table.hosts"
#' TODO: describe the table
#'
#' @param res an object of class \code{nosoiSim}
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with hosts table informations
#'
#' @export
#'
##
getTableHosts <- function(res, pop = "A") {
  return(getHostInfo(res, "table.hosts", pop))
}

#' @title Get table states
#'
#' @description
#' Get "table.state"
#' TODO: describe the table
#'
#' @param res an object of class \code{nosoiSim}
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with state table informations
#'
#' @export
#'
##
getTableState <- function(res, pop = "A") {
  return(getHostInfo(res, "table.state", pop))
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
  if (nosoiSim$popStructure == "none") return(NA)
  if (nosoiSim$popStructure == "discrete") return(df.meetTransmit[i, ]$current.in)
  if (nosoiSim$popStructure == "continuous") return(c(df.meetTransmit[i, ]$current.in.x, df.meetTransmit[i, ]$current.in.y))
  stop(paste0("Geographical structure ", nosoiSim$popStructure, " is not implemented."))
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
  if (nosoiSim$popStructure == "none") return(FALSE)
  if (nosoiSim$popStructure == "discrete") return(TRUE)
  if (nosoiSim$popStructure == "continuous") return(TRUE)
  stop(paste0("Geographical structure \"", nosoiSim$popStructure, "\" is not implemented."))
}

#' @title Param concatenator
#'
#' @description
#' Creates a \code{ParamHost} object.
#'
#' @inheritParams singleContinuous
#'
#' @return An object of class \code{ParamHost}
#'
#' @keywords internal
##

paramConstructor <- function(param.pExit, param.pMove, param.nContact, param.pTrans,
                             param.sdMove) {

  checkna <- function(param){return((length(param) == 1) && is.na(param))}

  if (checkna(param.pExit)) param.pExit <- NULL
  if (checkna(param.pMove)) param.pMove <- NULL
  if (checkna(param.nContact)) param.nContact <- NULL
  if (checkna(param.sdMove)) param.sdMove <- NULL
  if (checkna(param.pTrans)) param.pTrans <- NULL

  merged <- c(param.pTrans, param.pMove, param.nContact, param.pExit, param.sdMove)

  ParamHost <- merged[!duplicated(merged)]

  if(!is.null(ParamHost)) class(ParamHost) <- "ParamHost"

  return(ParamHost)
}
