#' @title Meet Transmit function
#'
#' @description
#' Perform the tasks related to meeting and transmission events.
#'
#' @param res an object of class \code{nosoiSim}.
#' @inheritParams newLine
#' @inheritParams newLineState
#'
#' @return The modified object res
#'
#' @keywords internal
##
meetTransmit <- function(res,
                         pres.time,
                         positions,
                         timeContactParsed, pTransParsed,
                         prefix.host, param.pExit, param.pMove, param.timeContact, param.pTrans,
                         param.moveDist) {

  keepHistory <- keepState(res)

  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

  df.meetTransmit <- res$table.hosts[active.hosts, c("hosts.ID", positions), with = FALSE]
  df.meetTransmit[, active.hosts:=hosts.ID]

  fun <- function(z) {
    timeContactParsed$vect(prestime = pres.time, z[, timeContactParsed$vectArgs, with = FALSE])
  }
  df.meetTransmit$number.contacts <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

  haveContact <- df.meetTransmit[["number.contacts"]] > 0
  df.meetTransmit <- df.meetTransmit[haveContact]
  active.hosts[active.hosts] <- haveContact # Update active hosts

  if (nrow(df.meetTransmit) > 0) {

    fun <- function(z) {
      pTransParsed$vect(prestime = pres.time, z[, pTransParsed$vectArgs, with = FALSE])
    }

    df.meetTransmit[, "Ptransmit"] <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"] #adds transmission probability to events
    df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Ptransmit"]] > 0] #discards event with probability 0

    if (nrow(df.meetTransmit) > 0) {

      df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

      df.meetTransmit[,"Trans"] <- drawBernouilli(df.meetTransmit[["Ptransmit"]]) #Draws K bernouillis with various probability (see function for more detail)

      df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Trans"]]] #Discards events with no realisation

      if (nrow(df.meetTransmit) >0) {
        table.temp <- vector("list", nrow(df.meetTransmit))
        if (keepHistory) table.state.temp <- vector("list", nrow(df.meetTransmit))
        for (i in 1:nrow(df.meetTransmit)) {

          res$N.infected <- res$N.infected + 1
          hosts.ID <- as.character(paste(prefix.host,res$N.infected,sep="-"))

          table.temp[[i]] <- newLine(hosts.ID,
                                     as.character(df.meetTransmit[i,]$active.hosts),
                                     infected.in = getPositionInfected(res, df.meetTransmit, i),
                                     pres.time,
                                     param.pExit, param.pMove, param.timeContact, param.pTrans,
                                     param.moveDist,
                                     current.environmental.value = df.meetTransmit[i,]$current.env.value)
          if (keepHistory) {
            table.state.temp[[i]] <- newLineState(hosts.ID,
                                                  state.pres = getPositionInfected(res, df.meetTransmit, i),
                                                  pres.time,
                                                  current.environmental.value = df.meetTransmit[i,]$current.env.value)
          }
        }

        res$table.hosts <- data.table::rbindlist(c(list(res$table.hosts),table.temp))
        data.table::setkey(res$table.hosts,hosts.ID)
        if (keepHistory) {
          res$table.state <- data.table::rbindlist(c(list(res$table.state),table.state.temp))
          data.table::setkey(res$table.state, "hosts.ID")
        }
      }
    }
  }

  return(res)
}
