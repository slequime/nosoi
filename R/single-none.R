#' Single-host without spatial features
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param timeContact function that gives the number of potential transmission events per unit of time.
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param pExit probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @details This function runs a single-host, without any spatial features, epidemiological simulation on a single core.
#' @details Th simulation stops either at the end of given time (specified by length) or when the number of hosts infected threshold is passed.
#' @export single_none
#' @import data.table

single_none <- function(length.sim,
                        max.infected,
                        init.individuals,
                        timeContact,
                        pTrans,
                        pExit,
                        prefix.host="H",
                        ...){

#Sanity check-------------------------------------------------------------
  if (is.na(length.sim) | length.sim <= 1) message("You must specify a length (in time units) for your simulation.")
  if (is.na(max.infected) | max.infected <= 1) message("You must specify a maximum number of infected hosts.")

  pTrans <- match.fun(pTrans)
  timeContact <- match.fun(timeContact)

  pTrans.param <- formalArgs(pTrans)[-1]
  n.pTrans.param <- length(pTrans.param)

  message("Starting the simulation")
  packageStartupMessage("Initializing ...", appendLF = FALSE)

#Creation of initial data ----------------------------------------------------------

  table.hosts <- iniTable(init.individuals,n.pTrans.param,pTrans.param,prefix.host)
  Host.count <- init.individuals

# Running the simulation ----------------------------------------
  packageStartupMessage(" running ...")
  pb <- txtProgressBar(min = 0, max = length.sim, style = 3,width=50)

for (pres.time in 1:length.sim){

  #Step 0: Active hosts ----------------------------------------------------------

  active.hosts = subset(table.hosts,active==1)$hosts.ID #active hosts

  if(length(active.hosts) > 0){
    exiting <- sample(c(TRUE,FALSE),length(active.hosts),replace=TRUE,prob=c(pExit,1-pExit))
    IDs = active.hosts[exiting]
    table.hosts[IDs, `:=` (out.time = as.numeric(pres.time),
                         active = 0)]

    active.hosts = subset(table.hosts, active==1)$hosts.ID #active hosts
  }

  if(length(active.hosts) == 0){break}

  #Step 1: Meeting & transmission ----------------------------------------------------

  for (j in active.hosts){
  number.contacts <- timeContact(1)

  for(contacts in 1:number.contacts)
    #is transmission occuring?

  x <- paste0("t=",pres.time-table.hosts[j, "inf.time"])

  x2 <- NULL

  for (i in pTrans.param){
    x1 <- paste0(i,"=",as.numeric(table.hosts[j, ..i]))
    x2 <- c(x2,x1)
  }

  option.trans <- paste0(c(x,x2),collapse = ",")

 Ptransmit <- eval(parse(text = paste("pTrans(", option.trans, ")")))

 transmit<-sample(c(TRUE,FALSE),1,prob=c(Ptransmit, 1-Ptransmit))

    if(transmit==TRUE){ #new host gets infected

      Host.count = Host.count+1
      hosts.ID = as.character(paste(prefix.host,Host.count,sep="-"))


      table.temp <- data.frame(matrix(0, ncol = (5+n.pTrans.param), nrow = 1))
      colnames(table.temp)=c("hosts.ID","inf.by","inf.time","out.time","active",pTrans.param)

      table.temp[1,"hosts.ID"] <- as.character(paste(prefix.host,Host.count,sep="-"))
      table.temp[1,"inf.by"] <- j
      table.temp[1,"inf.time"] <- pres.time
      table.temp[1,"out.time"] <- NA
      table.temp[1,"active"] <- 1
      for (i in pTrans.param){ table.temp[1,i] <- eval(parse(text = paste0(i,"(",1,")"))) }


      table.temp <- data.table::as.data.table(table.temp)
 data.table::setkey(table.temp, hosts.ID)


 table.hosts <- data.table::rbindlist(list(table.hosts,table.temp))
 data.table::setkey(table.hosts,hosts.ID)

  }
setTxtProgressBar(pb, pres.time)}
  if(Host.count > max.infected){break}
}
  packageStartupMessage(" done.")
  packageStartupMessage("The simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")

return(table.hosts)

}
