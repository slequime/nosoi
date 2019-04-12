#' Generates initial table to start the simulation (internal fonction)
#'
#' @description This function creates the initial table for the host.
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param init.structure State of the initially infected individuals.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param param.pExit list of exit probability function(s).
#' @param param.pMove list of movement probability function(s).
#' @param param.timeContact list of contact probability function(s).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @keywords internal

iniTable <- function(init.individuals,init.structure,prefix.host,param.pExit,param.pMove,param.timeContact,param.pTrans){

  list.init <- vector("list", init.individuals)

  for (indiv in 1:init.individuals) {
      list.init[[indiv]] <- newLine(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                    infected.by = paste(NA,indiv,sep="-"),
                                    infected.in = init.structure,
                                    time.is = 0,
                                    param.pExit = param.pExit,
                                    param.pMove = param.pMove,
                                    param.timeContact = param.timeContact,
                                    param.pTrans = param.pTrans)
  }

  table.hosts <- data.table::rbindlist(list.init)
  data.table::setkey(table.hosts, "hosts.ID")

  return(table.hosts)
}
