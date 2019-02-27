#' Generates initial table to start the simulation (internal fonction)
#'
#' @description This function creates the initial table for the host, with 5+number of parameters of the transmission probability function paramters, and init.individuals row(s).
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param n.pExit.param number of parameters for exit probability computation (impacts number of columns).
#' @param param.pExit list of exit probability function(s).
#' @param n.pTrans.param number of parameters for transmission probability computation (impacts number of columns).
#' @param param.pTrans list of transmission probability function(s).
#'
#' @keywords internal

iniTable <- function(init.individuals,prefix.host, n.pExit.param,param.pExit,n.pTrans.param, param.pTrans){
  #Creation of initial data ----------------------------------------------------------

  # table.hosts <- data.frame(matrix(0, ncol = (5+n.pTrans.param+n.pExit.param), nrow = init.individuals))
  # colnames(table.hosts)<-c("hosts.ID","inf.by","inf.time","out.time","active",names(param.pExit),names(param.pTrans))
  list.init <- vector("list", init.individuals)

  for(indiv in 1:init.individuals){

    list.init[[indiv]] <- newLine(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                  infected.by = paste(NA,indiv,sep="-"),
                                  time.is = 0,
                                  n.pExit.param = n.pExit.param,
                                  param.pExit = param.pExit,
                                  n.pTrans.param = n.pTrans.param,
                                  param.pTrans = param.pTrans)

    # table.hosts[indiv,"hosts.ID"] <- paste(prefix.host,indiv,sep="-")
    # table.hosts[indiv,"inf.by"] <- paste(NA,indiv,sep="-")
    # table.hosts[indiv,"inf.time"] <- 0
    # table.hosts[indiv,"out.time"] <- NA
    # table.hosts[indiv,"active"] <- 1
    #
    # Sampled.param.pTrans <- sapply(param.pTrans, function(x) sapply(1, x))
    #
    # for (i in names(param.pTrans)){ table.hosts[indiv,i] = Sampled.param.pTrans[i] }
    # if(n.pExit.param > 0){
    #   Sampled.param.pExit <- sapply(param.pExit, function(x) sapply(1, x))
    #   for (j in names(param.pExit)){ table.hosts[indiv,j] = Sampled.param.pExit[j] }
    # }
  }
  table.hosts <- data.table::rbindlist(list.init)

  # table.hosts <- data.table::as.data.table(table.hosts)
  data.table::setkey(table.hosts, "hosts.ID")
  return(table.hosts)
}
