#' Single-host without structured host population
#'
#' @description This function runs a single-host transmission chain simulation, without any spatial features. The simulation stops either at
#' the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected) is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param timeContact function that gives the number of potential transmission events per unit of time.
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pExit.
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @export singleNone

singleNone <- function(length.sim,
                       max.infected,
                       init.individuals,
                       timeContact,
                       pTrans,
                       param.pTrans,
                       pExit,
                       param.pExit,
                       prefix.host="H",
                       progress.bar=TRUE,
                       ...){

#Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  if (is.na(length.sim) | length.sim <= 1) stop("You must specify a length (in time units) for your simulation (bigger than 1).")
  if (is.na(max.infected) | max.infected <= 1) stop("You must specify a maximum number of infected hosts (bigger than 1).")
  if (is.na(init.individuals) | init.individuals < 1 | !init.individuals%%1==0) stop("The transmission chain should be started by 1 or more (integer) individuals")

  if (! is.function(timeContact)) stop("Contact probability should be a function of time.")
  if (! is.function(pTrans)) stop("Transmission probability should be a function of time.")
  if (! is.function(pExit)) stop("Exit probability should be a function of time.")

  timeContact <- match.fun(timeContact)

  pTrans <- match.fun(pTrans)

  if (is.list(param.pTrans)){
    pTrans.param <- formalArgs(pTrans)[-1]
    if(! all(names(param.pTrans) %in% pTrans.param)) stop("Parameter name in param.pTrans should match the name used in pTrans.")
    n.pTrans.param <- length(pTrans.param)
  }

  pExit <- match.fun(pExit)

  if(length(formalArgs(pExit)) == 1){pExit.type <- "simple"
  n.pExit.param <- 0}
  if(length(formalArgs(pExit)) > 1){pExit.type <- "complex"

  if (! is.na(param.pExit)){
    pExit.param <- formalArgs(pExit)[-1]
    if(! all(names(param.pExit) %in% pExit.param)) stop("Parameter name in param.pExit should match the name used in pExit.")
    n.pExit.param <- length(pExit.param)
    pExit.constant <- FALSE
    }else(stop("There is a probleme with your exit probability function: you should provide a parameter list."))
  }

  message("Starting the simulation")
  packageStartupMessage("Initializing ...", appendLF = FALSE)

#Creation of initial data ----------------------------------------------------------

  table.hosts <- iniTable(init.individuals,prefix.host,n.pExit.param,param.pExit, n.pTrans.param, param.pTrans,...)
  Host.count <- init.individuals

# Running the simulation ----------------------------------------
  packageStartupMessage(" running ...")
if( progress.bar==TRUE){pb <- txtProgressBar(min = 0, max = length.sim, style = 3, width=50)}

for (pres.time in 1:length.sim){

  #Step 0: Active hosts ----------------------------------------------------------

  active.hosts = subset(table.hosts,active==1)$hosts.ID #active hosts
  if(length(active.hosts) > 0){

    if(pExit.type == "simple"){
      p.exit.values = pExit(pres.time-table.hosts[active.hosts]$inf.time)
      exiting = drawBernouilli(p.exit.values)
    }

    if(pExit.type == "complex"){
      p.exit.values = NULL
      for(j in active.hosts){
        x2 <- NULL
        for (i in names(param.pExit)){x2 <- c(x2,paste0(i,"=",as.numeric(table.hosts[j, ..i])))}

        option.trans <- paste0(c(paste0("t=",pres.time-table.hosts[j, "inf.time"]),x2),collapse = ",")
        p.exit.values = c(p.exit.values,eval(parse(text = paste("pExit(", option.trans, ")"))))
      }

      exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
    }
  }

  IDs = active.hosts[exiting]
  table.hosts[IDs, `:=` (out.time = as.numeric(pres.time),
                         active = 0)]

  active.hosts = subset(table.hosts, active==1)$hosts.ID #active hosts
  if(length(active.hosts) == 0){break}
  #Step 1: Meeting & transmission ----------------------------------------------------

  number.contacts <- timeContact(length(active.hosts))

  df.meetTransmit = data.table(active.hosts,number.contacts)
  data.table::setkey(df.meetTransmit,active.hosts)

  df.meetTransmit = subset(df.meetTransmit, number.contacts > 0)

  if(nrow(df.meetTransmit) >0){

  results.Ptransmit=NULL
    for(j in df.meetTransmit$active.hosts){

      x2 <- NULL
      for (i in names(param.pTrans)){x2 <- c(x2,paste0(i,"=",as.numeric(table.hosts[j, ..i])))}

      option.trans <- paste0(c(paste0("t=",pres.time-table.hosts[j, "inf.time"]),x2),collapse = ",")
      results.Ptransmit = c(results.Ptransmit,eval(parse(text = paste("pTrans(", option.trans, ")"))))
    }

  df.meetTransmit$Ptransmit = results.Ptransmit #adds transmission probability to events
  df.meetTransmit <- subset(df.meetTransmit, Ptransmit > 0) #discards event with probability 0

  if(nrow(df.meetTransmit) >0){

  df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

  df.meetTransmit$Trans <- drawBernouilli(df.meetTransmit$Ptransmit) #Draws K bernouillis with various probability (see function for more detail)

  df.meetTransmit = subset(df.meetTransmit, Trans == TRUE) #Discards events with no realisation

  if(nrow(df.meetTransmit) >0){
  for(i in 1:nrow(df.meetTransmit)){

    Host.count <- Host.count+1
    hosts.ID <- as.character(paste(prefix.host,Host.count,sep="-"))

    table.temp <- newLine(hosts.ID,as.character(df.meetTransmit[i,]$active.hosts),pres.time,n.pExit.param,param.pExit,n.pTrans.param,param.pTrans)

    table.hosts <- data.table::rbindlist(list(table.hosts,table.temp))
    data.table::setkey(table.hosts,hosts.ID)
  }}}}


if( progress.bar==TRUE){setTxtProgressBar(pb, pres.time)}
  if(Host.count > max.infected){break}
}
  packageStartupMessage(" done.")
  packageStartupMessage("The simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")

return(table.hosts)

}
