#' @title Meet Transmit function
#'
#' @description
#' Perform the tasks related to meeting and transmission events.
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @inheritParams newLine
#' @inheritParams newLineState
#'
#' @return A table with active transmission events
#'
#' @keywords internal
##
meetTransmit <- function(res,
                         pres.time,
                         positions,
                         nContactParsed, pTransParsed) {
  #To avoids notes (use of dplyr functions)
  hosts.ID <- NULL

  active.hosts <- res$table.hosts[["active"]] #active hosts (boolean vector)

  df.meetTransmit <- res$table.hosts[active.hosts, c("hosts.ID", positions), with = FALSE]
  df.meetTransmit[, active.hosts:=hosts.ID]

  df.meetTransmit$number.contacts <- applyFunctionToHosts(res, pres.time, nContactParsed, active.hosts)

  haveContact <- df.meetTransmit[["number.contacts"]] > 0
  df.meetTransmit <- df.meetTransmit[haveContact]
  active.hosts[active.hosts] <- haveContact # Update active hosts

  if (nrow(df.meetTransmit) > 0) {

    df.meetTransmit[, "Ptransmit"] <- applyFunctionToHosts(res, pres.time, pTransParsed, active.hosts) #adds transmission probability to events
    df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Ptransmit"]] > 0] #discards event with probability 0

    if (nrow(df.meetTransmit) > 0) {
      df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

      df.meetTransmit[,"Trans"] <- drawBernouilli(df.meetTransmit[["Ptransmit"]]) #Draws K bernouillis with various probability (see function for more detail)

      df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Trans"]]] #Discards events with no realisation
    }

  }
  return(df.meetTransmit)
}

#' @title Write newly infected function
#'
#' @description
#' Perform the tasks related to creating new lines of hosts in relevant tables.
#'
#' @param df.meetTransmit table with active transmission events, coming from \code{meetTransmit}
#' @param res an object of class \code{nosoiSimOne}.
#' @inheritParams newLine
#' @inheritParams newLineState
#'
#' @return The modified object res
#'
#' @keywords internal
##

writeInfected <- function(df.meetTransmit, res,
                          pres.time, ParamHost) {

  keepHistory <- keepState(res)

  if (nrow(df.meetTransmit) > 0) {

    table.temp <- vector("list", nrow(df.meetTransmit))
    if (keepHistory) table.state.temp <- vector("list", nrow(df.meetTransmit))
    for (i in 1:nrow(df.meetTransmit)) {

      res$N.infected <- res$N.infected + 1
      hosts.ID <- as.character(paste(res$prefix.host, res$N.infected, sep="-"))

      table.temp[[i]] <- newLine(hosts.ID,
                                 as.character(df.meetTransmit[[i, "active.hosts"]]),
                                 infected.in = getPositionInfected(res, df.meetTransmit, i),
                                 pres.time,
                                 ParamHost,
                                 current.environmental.value = df.meetTransmit[[i, "current.env.value"]],
                                 current.cell.number.raster = df.meetTransmit[[i, "current.cell.raster"]],
                                 current.count.A = df.meetTransmit[[i, "host.count.A"]],
                                 current.count.B = df.meetTransmit[[i, "host.count.B"]])
      if (keepHistory) {
        table.state.temp[[i]] <- newLineState(hosts.ID,
                                              state.pres = getPositionInfected(res, df.meetTransmit, i),
                                              pres.time,
                                              current.environmental.value = df.meetTransmit[[i, "current.env.value"]],
                                              current.cell.number.raster = df.meetTransmit[[i, "current.cell.raster"]])
      }
    }

    res$table.hosts <- data.table::rbindlist(c(list(res$table.hosts),table.temp))
    data.table::setkey(res$table.hosts,hosts.ID)
    if (keepHistory) {
      res$table.state <- data.table::rbindlist(c(list(res$table.state),table.state.temp))
      data.table::setkey(res$table.state, "hosts.ID")
    }

  }
  return(res)
}
