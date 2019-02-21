#' Generates initial table to start the simulation (internal fonction)
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param n.pTrans.param number of parameters for transmission probability computation (impacts number of columns).
#' @param pTrans.param transmission probability function parameters' name.
#'
#' @details This function creates the initial table for the host, with 5+number of paramters of the transmission probability function paramters, and init.individuals row(s).
#' @export iniTable
#' @import data.table

iniTable <- function(init.individuals,n.pTrans.param,pTrans.param){

#Creation of initial data ----------------------------------------------------------
table.hosts <- data.frame(matrix(0, ncol = (5+n.pTrans.param), nrow = init.individuals))
colnames(table.hosts)<-c("hosts.ID","inf.by","inf.time","out.time","active",pTrans.param)

for(indiv in 1:init.individuals){

  table.hosts[indiv,"hosts.ID"] <- paste("H",indiv,sep="-")
  table.hosts[indiv,"inf.by"] <- paste("unkown")
  table.hosts[indiv,"inf.time"] <- 0
  table.hosts[indiv,"out.time"] <- NA
  table.hosts[indiv,"active"] <- 1
  for (i in pTrans.param){ table.hosts[indiv,i] = eval(parse(text = paste0(i,"(",1,")"))) }
}

table.hosts <- data.table::as.data.table(table.hosts)
data.table::setkey(table.hosts, hosts.ID)
return(table.hosts)
}
