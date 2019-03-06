#' @title Creates a new line to be added to the table when new host is infected (internal fonction)
#'
#' @description
#' This function creates a new line for the table,
#' with 5+number of parameters of the transmission probability function paramters,
#' and init.individuals row(s).
#' The lines are to be bounded with \code{\link{data.table::rbindlist}}.
#'
#' @param hosts.ID unique ID for the new host
#' @param infected.by unique ID of host that transmits to the new one
#' @param time.is time in the simulation, when the infection takes place
#' @param param.pExit list of exit probability function(s).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @return a list with the new line to add.
#'
#' @keywords internal

newLine <- function(hosts.ID,infected.by,infected.in,time.is,param.pExit,param.pMove,param.pTrans) {
  if (is.na(param.pExit)) param.pExit <- NULL
  if (is.na(param.pMove)) param.pMove <- NULL
  if (is.na(infected.in)) infected.in <- NULL

  return(c(hosts.ID = hosts.ID,
           inf.by = infected.by,
           inf.in = infected.in,
           inf.time = time.is,
           out.time = NA_real_,
           active = 1,
           as.list(sapply(param.pTrans, function(x) x(1))),
           as.list(sapply(param.pMove, function(x) x(1))),
           as.list(sapply(param.pExit, function(x) x(1)))
  )
  )
}
