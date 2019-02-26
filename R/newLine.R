#' Creates a new line to be added to the table when new host is infected (internal fonction)
#'
#' @description This function creates the initial table for the host, with 5+number of parameters of the transmission probability function paramters, and init.individuals row(s).
#'
#' @param hosts.ID unique ID for the new host
#' @param infected.by unique ID of host that transmits to the new one
#' @param time.is time in the simulation, when the infection takes place
#' @param n.pExit.param number of parameters for exit probability computation (impacts number of columns).
#' @param param.pExit list of exit probability function(s).
#' @param n.pTrans.param number of parameters for transmission probability computation (impacts number of columns).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @keywords internal

newLine <- function(hosts.ID,infected.by,time.is,n.pExit.param,param.pExit,n.pTrans.param,param.pTrans) {

  table.temp <- data.frame(matrix(0, ncol = (5+n.pTrans.param+n.pExit.param), nrow = 1))
  colnames(table.temp)<-c("hosts.ID","inf.by","inf.time","out.time","active",names(param.pExit),names(param.pTrans))

  table.temp[1,"hosts.ID"] <- hosts.ID
  table.temp[1,"inf.by"] <- infected.by
  table.temp[1,"inf.time"] <- time.is
  table.temp[1,"out.time"] <- NA
  table.temp[1,"active"] <- 1

  Sampled.param <- sapply(param.pTrans, function(x) sapply(1, x))

  for (i in names(param.pTrans)) { table.temp[1,i] <- Sampled.param[i] }

  if (n.pExit.param > 0) {
    Sampled.param.pExit <- sapply(param.pExit, function(x) sapply(1, x))
    for (j in names(param.pExit)) { table.temp[1,j] = Sampled.param.pExit[j] }
  }

  table.temp <- data.table::as.data.table(table.temp)

  return(table.temp)
}



