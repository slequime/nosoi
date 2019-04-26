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

#' @title Output wrapper
#'
#' @description
#' Wraps up the nosoi output
#'
#' @param Host.count number of infected hosts
#' @param pres.time current time of the simulation
#' @param table.hosts data.table of hosts
#' @param state.archive data.table of hosts movement
#'
#'
#' @keywords internal
##

outputWrapper <- function(Host.count, pres.time, table.hosts, state.archive) {

return(list(total.time = pres.time,
            N.infected = Host.count,
            table.hosts = table.hosts,
            table.state = state.archive))
}
