#' Single-host, no spatial features and on a single core
#'
#' @param length specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals if simulation should be run in parallel or not.
#' @param time_contact function that gives the number of potential transmission events per unit of time.
#' @param p_trans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param p_exit probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @details This function runs a single-host, without any spatial features, epidemiological simulation on a single core.
#' @details Th simulation stops either at the end of given time (specified by length) or when the number of hosts infected threshold is passed.
#' @export single_none_parallel
#' @import data.table

single_none_parallel <- function(length,
                                 max.infected,
                                 init.individuals,
                                 time_contact,
                                 p_trans,
                                 p_exit,
                                 ...){

#Function check-------------------------------------------------------------

  p_trans <- match.fun(p_trans)
  time_contact <- match.fun(time_contact)

  p_trans.param <- formalArgs(p_trans)[-1]
  n.p_trans.param <- length(p_trans.param)

  message("Starting the simulation")
  packageStartupMessage("Initializing ...", appendLF = FALSE)

#Creation of initial data ----------------------------------------------------------
  table.hosts <- data.frame(matrix(0, ncol = (5+n.p_trans.param), nrow = init.individuals))
  colnames(table.hosts)<-c("hosts.ID","inf.by","inf.time","out.time","active",p_trans.param)

  for(indiv in 1:init.individuals){

    table.hosts[indiv,"hosts.ID"] <- paste("H",indiv,sep="-")
    table.hosts[indiv,"inf.by"] <- paste("unkown")
    table.hosts[indiv,"inf.time"] <- 0
    table.hosts[indiv,"out.time"] <- NA
    table.hosts[indiv,"active"] <- 1
  for (i in p_trans.param){ table.hosts[indiv,i] = eval(parse(text = paste0(i,"(",1,")"))) }
  }

  table.hosts <- data.table::as.data.table(table.hosts)
  data.table::setkey(table.hosts, hosts.ID)
  Host.count <- init.individuals

# Running the simulation ----------------------------------------
  packageStartupMessage(" running ...")
  pb <- txtProgressBar(min = 0, max = length, style = 3,width=50)

for (pres.time in 1:length){

  #Step 0: Active hosts ----------------------------------------------------------

  active.hosts = subset(table.hosts,active==1)$hosts.ID #active hosts

  if(length(active.hosts) > 0){
    exiting <- sample(c(TRUE,FALSE),length(active.hosts),replace=TRUE,prob=c(p_exit,1-p_exit))
    IDs = active.hosts[exiting]
    table.hosts[IDs, `:=` (out.time = as.numeric(pres.time),
                         active = 0)]

    active.hosts = subset(table.hosts, active==1)$hosts.ID #active hosts
  }

  if(length(active.hosts) == 0){break}

  #Step 1: Meeting & transmission ----------------------------------------------------

  for (j in active.hosts){
  number.contacts <- time_contact(1)

  for(contacts in 1:number.contacts)
    #is transmission occuring?

  x <- paste0("t=",pres.time-table.hosts[j, "inf.time"])

  x2 <- NULL

  for (i in p_trans.param){
    x1 <- paste0(i,"=",as.numeric(table.hosts[j, ..i]))
    x2 <- c(x2,x1)
  }

  option.trans <- paste0(c(x,x2),collapse = ",")

 Ptransmit <- eval(parse(text = paste("p_trans(", option.trans, ")")))

 transmit<-sample(c(TRUE,FALSE),1,prob=c(Ptransmit, 1-Ptransmit))

    if(transmit==TRUE){ #new host gets infected

      Host.count = Host.count+1
      hosts.ID = as.character(paste("H",Host.count,sep="-"))


      table.temp <- data.frame(matrix(0, ncol = (5+n.p_trans.param), nrow = 1))
      colnames(table.temp)=c("hosts.ID","inf.by","inf.time","out.time","active",p_trans.param)

      table.temp[1,"hosts.ID"] <- as.character(paste("H",Host.count,sep="-"))
      table.temp[1,"inf.by"] <- j
      table.temp[1,"inf.time"] <- pres.time
      table.temp[1,"out.time"] <- NA
      table.temp[1,"active"] <- 1
      for (i in p_trans.param){ table.temp[1,i] <- eval(parse(text = paste0(i,"(",1,")"))) }


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
