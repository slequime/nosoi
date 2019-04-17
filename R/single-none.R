#' Single-host without structured host population
#'
#' @description This function runs a single-host transmission chain simulation, without any spatial features. The simulation stops either at
#' the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected) is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param timeContact function that gives the number of potential transmission events per unit of time.
#' @param param.timeContact parameter names (list of functions) for param.timeContact.
#' @param timeDep.timeContact is timeContact dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pExit.
#' @param timeDep.pTrans is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param timeDep.pExit is pExit dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param progress.bar if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @export singleNone

singleNone <- function(length.sim,
                       max.infected,
                       init.individuals,
                       timeContact,
                       param.timeContact,
                       timeDep.timeContact=FALSE,
                       pTrans,
                       param.pTrans,
                       timeDep.pTrans=FALSE,
                       pExit,
                       param.pExit,
                       timeDep.pExit=FALSE,
                       prefix.host="H",
                       progress.bar=TRUE,
                       print.step=10,
                       ...){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  if (is.na(length.sim) | length.sim <= 1) stop("You must specify a length (in time units) for your simulation (bigger than 1).")
  if (is.na(max.infected) | max.infected <= 1) stop("You must specify a maximum number of infected hosts (bigger than 1).")
  if (is.na(init.individuals) | init.individuals < 1 | !init.individuals%%1==0) stop("The transmission chain should be started by 1 or more (integer) individuals")

  if (! is.function(timeContact)) stop("Contact probability should be a function of time.")
  if (! is.function(pTrans)) stop("Transmission probability should be a function of time.")
  if (! is.function(pExit)) stop("Exit probability should be a function of time.")

  #Parsing timeContact
  timeContact <- match.fun(timeContact)

  if (timeDep.timeContact == FALSE) {
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=FALSE,timeDep=FALSE)
  }

  if (timeDep.timeContact == TRUE) {
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),'prestime'))==FALSE) stop("timeContact should have 'prestime' as a variable. timeDep.timeContact == TRUE.")
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=FALSE,timeDep=TRUE)
  }

  #Parsing pTrans

  if (timeDep.pTrans == FALSE) {
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=FALSE,timeDep=FALSE)
  }

  if (timeDep.pTrans == TRUE) {
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),'prestime'))==FALSE) stop("pTrans should have 'prestime' as a variable. timeDep.pTrans == TRUE.")
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=FALSE,timeDep=TRUE)
  }

  #Parsing pExit
  if (timeDep.pExit == FALSE) {
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=FALSE,timeDep=FALSE)
  }

  if (timeDep.pExit == TRUE) {
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),'prestime'))==FALSE) stop("pExit should have 'prestime' as a variable. timeDep.pExit == TRUE.")
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=FALSE,timeDep=TRUE)
  }

  # Init
  message("Starting the simulation")
  message("Initializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  table.hosts <- iniTable(init.individuals, NA, prefix.host, param.pExit, param.pMove = NA, param.timeContact, param.pTrans, param.moveDist=NA)
  Host.count <- init.individuals

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    active.hosts <- table.hosts[["active"]] == 1 #active hosts (boolean vector)
    if (any(active.hosts)) {

      fun <- function(z) {
        pExitParsed$vect(prestime = pres.time, z[, pExitParsed$vectArgs, with = FALSE])
      }
      p.exit.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

      exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
    }
    # }

    exiting.full <- active.hosts
    exiting.full[exiting.full] <- exiting

    table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                    active = 0)]

    active.hosts[active.hosts] <- !exiting # Update active hosts

    if (!any(active.hosts)) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------
    df.meetTransmit <- table.hosts[active.hosts, c("hosts.ID")]
    df.meetTransmit[, active.hosts:=hosts.ID]

    fun <- function(z) {
      timeContactParsed$vect(prestime = pres.time, z[, timeContactParsed$vectArgs, with = FALSE])
    }
    timeContact.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

    df.meetTransmit$number.contacts <- timeContact.values

    haveContact <- df.meetTransmit[["number.contacts"]] > 0
    df.meetTransmit <- df.meetTransmit[haveContact]
    active.hosts[active.hosts] <- haveContact # Update active hosts

    if (nrow(df.meetTransmit) > 0) {

      fun <- function(z) {
        pTransParsed$vect(prestime = pres.time, z[, pTransParsed$vectArgs, with = FALSE])
      }

      df.meetTransmit[, "Ptransmit"] <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"] #adds transmission probability to events
      df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Ptransmit"]] > 0] #discards event with probability 0

      if (nrow(df.meetTransmit) > 0) {

        df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

        df.meetTransmit[,"Trans"] <- drawBernouilli(df.meetTransmit[["Ptransmit"]]) #Draws K bernouillis with various probability (see function for more detail)

        df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Trans"]]] #Discards events with no realisation

        if (nrow(df.meetTransmit) >0) {
          table.temp <- vector("list", nrow(df.meetTransmit))
          for (i in 1:nrow(df.meetTransmit)) {

            Host.count <- Host.count+1
            hosts.ID <- as.character(paste(prefix.host,Host.count,sep="-"))

            table.temp[[i]] <- newLine(hosts.ID, as.character(df.meetTransmit[i,]$active.hosts), NA, pres.time, param.pExit, param.pMove=NA,param.timeContact, param.pTrans,param.moveDist=NA)
          }

          table.hosts <- data.table::rbindlist(c(list(table.hosts),table.temp))
          data.table::setkey(table.hosts,hosts.ID)
        }
      }
    }

    if (progress.bar == TRUE & pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count," (",round((Host.count/max.infected)*100,digits=0),"% of maximum infected hosts).")}
    if (Host.count > max.infected) {break}
  }
  message(" done.")
  message("The simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")

  nosoi.output <- list()

  nosoi.output[["total.time"]] <- pres.time
  nosoi.output[["N.infected"]] <- Host.count
  nosoi.output[["table.hosts"]] <- table.hosts

  return(nosoi.output)
}
