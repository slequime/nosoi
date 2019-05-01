#' @title Meet Transmit function
#'
#' @description
#' Perform the tasks related to meeting and transmission events.
#'
#' @param res an object of class \code{nosoiSim}.
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
                         timeContactParsed, pTransParsed) {

  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

  df.meetTransmit <- res$table.hosts[active.hosts, c("hosts.ID", positions), with = FALSE]
  df.meetTransmit[, active.hosts:=hosts.ID]

  df.meetTransmit$number.contacts <- applyFunctionToHosts(res, pres.time, timeContactParsed, active.hosts)

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
