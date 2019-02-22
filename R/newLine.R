#' Creates a new line to be added to the table when new host is infected (internal fonction)
#'
#' @param hosts.ID unique ID for the new host
#' @param infected.by unique ID of host that transmits to the new one
#' @param time.is time in the simulation, when the infection takes place
#' @param n.pTrans.param number of parameters for transmission probability computation (impacts number of columns).
#' @param pTrans.param transmission probability function parameters' name.
#'
#' @details This function creates the initial table for the host, with 5+number of parameters of the transmission probability function paramters, and init.individuals row(s).
#' @export newLine
#' @import data.table

newLine <- function(hosts.ID,infected.by,time.is,n.pTrans.param,param.pTrans){
  #new host gets infected
  table.temp <- data.frame(matrix(0, ncol = (5+n.pTrans.param), nrow = 1))
  colnames(table.temp) <- c("hosts.ID","inf.by","inf.time","out.time","active",names(param.pTrans))

  table.temp[1,"hosts.ID"] <- hosts.ID
  table.temp[1,"inf.by"] <- infected.by
  table.temp[1,"inf.time"] <- time.is
  table.temp[1,"out.time"] <- NA
  table.temp[1,"active"] <- 1

  Sampled.param <- sapply(param.pTrans, function(x) sapply(1, x))

  for (i in names(param.pTrans)){ table.temp[1,i] <- Sampled.param[i] }

  table.temp <- data.table::as.data.table(table.temp)
  data.table::setkey(table.temp, hosts.ID)

return(table.temp)
}



