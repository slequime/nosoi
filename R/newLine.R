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
#' @param n.pExit.param number of parameters for exit probability computation (impacts number of columns).
#' @param param.pExit list of exit probability function(s).
#' @param n.pTrans.param number of parameters for transmission probability computation (impacts number of columns).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @return a list with the new line to add.
#'
#' @keywords internal

newLine <- function(hosts.ID,infected.by,time.is,n.pExit.param,param.pExit,n.pTrans.param,param.pTrans) {
  if (is.na(param.pExit)) param.pExit <- NULL
  return(c(hosts.ID = hosts.ID,
           inf.by = infected.by,
           inf.time = time.is,
           out.time = NA_real_,
           active = 1,
           as.list(sapply(param.pTrans, function(x) x(1))),
           as.list(sapply(param.pExit, function(x) x(1)))
  )
  )

  # table.temp <- data.frame(matrix(0, ncol = (5+n.pTrans.param+n.pExit.param), nrow = 1))
  # colnames(table.temp)<-c("hosts.ID","inf.by","inf.time","out.time","active",names(param.pExit),names(param.pTrans))
  #
  # table.temp[1,"hosts.ID"] <- hosts.ID
  # table.temp[1,"inf.by"] <- infected.by
  # table.temp[1,"inf.time"] <- time.is
  # table.temp[1,"out.time"] <- NA
  # table.temp[1,"active"] <- 1
  #
  # Sampled.param <- sapply(param.pTrans, function(x) sapply(1, x))
  #
  # for (i in names(param.pTrans)) { table.temp[1,i] <- Sampled.param[i] }
  #
  # if (n.pExit.param > 0) {
  #   Sampled.param.pExit <- sapply(param.pExit, function(x) sapply(1, x))
  #   for (j in names(param.pExit)) { table.temp[1,j] = Sampled.param.pExit[j] }
  # }
  #
  # table.temp <- data.table::as.data.table(table.temp)
  #
  # return(table.temp)
}



