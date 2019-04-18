#' @title Progress bar
#'
#' @description
#' Echos the state of the simulation at any given time step provided by the user.
#'
#' @param Host.count
#' @param pres.time
#' @param print.step
#' @param length.sim
#' @param max.infected
#'
#'
#' @seealso \code{\link{nosoiSim}}
##

progressMessage <- function(Host.count, pres.time, print.step, length.sim, max.infected) {
  if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count," (",round((Host.count/max.infected)*100,digits=0),"% of maximum infected hosts).")}
}

#' @title End message
#'
#' @description
#' Message that ends the simulation
#'
#' @param Host.count
#' @param pres.time
#'
#'
#' @seealso \code{\link{nosoiSim}}
##

endMessage <- function(Host.count, pres.time) {
  message("done. \nThe simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")
  }

#' @title Output wrapper
#'
#' @description
#' Wraps up the nosoi output
#'
#' @param Host.count
#' @param pres.time
#' @param table.hosts
#' @param state.archive
#'
#'
#' @seealso \code{\link{nosoiSim}}
##

outputWrapper <- function(Host.count, pres.time, table.hosts, state.archive) {

return(list(total.time = pres.time,
            N.infected = Host.count,
            table.hosts = table.hosts,
            table.state = state.archive))
}
