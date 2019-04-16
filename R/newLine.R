#' @title Creates a new line to be added to the table when new host is infected (internal fonction)
#'
#' @description
#' This function creates a new line for the table.
#' The lines are to be bounded with \code{\link[data.table]{rbindlist}}.
#'
#' @param hosts.ID unique ID for the new host
#' @param infected.by unique ID of host that transmits to the new one
#' @param infected.in state in which the host was infected
#' @param current.in state in which the host currently is
#' @param time.is time in the simulation, when the infection takes place
#' @param param.pExit list of exit probability function(s).
#' @param param.pMove list of movement probability function(s).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @return a list with the new line to add.
#'
#' @keywords internal

newLine <- function(hosts.ID,infected.by,infected.in,time.is,param.pExit,param.pMove,param.timeContact,param.pTrans,current.environmental.value=NA) {
  if (is.na(param.pExit)) param.pExit <- NULL
  if (is.na(param.pMove)) param.pMove <- NULL
  if (is.na(param.timeContact)) param.timeContact <- NULL
  if (is.na(current.environmental.value)) current.environmental.value <- NULL

  if (length(infected.in) == 1){
    if (is.na(infected.in)) infected.in <- NULL
  return(c(hosts.ID = hosts.ID,
           inf.by = infected.by,
           inf.in = infected.in,
           current.in = infected.in,
           inf.time = time.is,
           out.time = NA_real_,
           active = 1,
           as.list(sapply(param.pTrans, function(x) x(1))),
           as.list(sapply(param.pMove, function(x) x(1))),
           as.list(sapply(param.timeContact, function(x) x(1))),
           as.list(sapply(param.pExit, function(x) x(1)))
  )
  )}

  if (length(infected.in) == 2){

    return(c(hosts.ID = hosts.ID,
             inf.by = infected.by,
             inf.in.x = infected.in[1],
             inf.in.y = infected.in[2],
             current.in.x = infected.in[1],
             current.in.y = infected.in[2],
             current.env.value = current.environmental.value,
             inf.time = time.is,
             out.time = NA_real_,
             active = 1,
             as.list(sapply(param.pTrans, function(x) x(1))),
             as.list(sapply(param.pMove, function(x) x(1))),
             as.list(sapply(param.timeContact, function(x) x(1))),
             as.list(sapply(param.pExit, function(x) x(1)))
    )
    )}
}
