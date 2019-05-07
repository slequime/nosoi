#' @title Progress bar
#'
#' @description
#' Echos the state of the simulation at any given time step provided by the user.
#'
#' @param Host.count.A number of infected hosts of host A.
#' @param Host.count.B number of infected hosts of host B.
#' @param pres.time current time of the simulation
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
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

nosoiSimConstructor <- function(N.infected, total.time, table.hosts, table.state, prefix.host,
                                type = c("singleNone", "singleDiscrete", "singleContinuous","dualNone","dualDiscrete", "dualContinuous")) {

  type <- match.arg(type)

  res <- list(total.time = total.time,
              N.infected = N.infected,
              table.hosts = table.hosts,
              table.state = table.state,
              prefix.host = prefix.host,
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
  if (nosoiSim$type == "dualNone") return(NA)
  if (nosoiSim$type == "dualDiscrete") return(df.meetTransmit[i, ]$current.in)
  if (nosoiSim$type == "dualContinuous") return(c(df.meetTransmit[i, ]$current.in.x, df.meetTransmit[i, ]$current.in.y))
  stop(paste("This type",nosoiSim$type,"is not implemented yet"))
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
  if (nosoiSim$type == "dualNone") return(FALSE)
  if (nosoiSim$type == "dualDiscrete") return(TRUE)
  if (nosoiSim$type == "dualContinuous") return(TRUE)
  stop(paste("This type",nosoiSim$type,"is not implemented yet"))
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

  class(ParamHost) <- "ParamHost"

  return(ParamHost)
}
