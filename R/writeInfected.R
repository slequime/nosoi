#' @title Write newly infected function
#'
#' @description
#' Perform the tasks related to creating new lines of hosts in relevant tables.
#'
#' @param df.meetTransmit table with active transmission events, coming from \code{meetTransmit}
#' @param res an object of class \code{nosoiSim}.
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
                                   as.character(df.meetTransmit[i,]$active.hosts),
                                   infected.in = getPositionInfected(res, df.meetTransmit, i),
                                   pres.time,
                                   ParamHost,
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
  return(res)
}
