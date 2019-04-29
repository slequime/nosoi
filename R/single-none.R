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
                       pExit,
                       param.pExit,
                       timeDep.pExit=FALSE,
                       timeContact,
                       param.timeContact,
                       timeDep.timeContact=FALSE,
                       pTrans,
                       param.pTrans,
                       timeDep.pTrans=FALSE,
                       prefix.host="H",
                       progress.bar=TRUE,
                       print.step=10,
                       ...){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  # if (! is.function(timeContact)) stop("Contact probability should be a function of time.")
  # if (! is.function(pTrans)) stop("Transmission probability should be a function of time.")
  # if (! is.function(pExit)) stop("Exit probability should be a function of time.")

  #Parsing timeContact
  timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),timeDep=timeDep.timeContact)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),timeDep=timeDep.pTrans)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),timeDep=timeDep.pExit)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(N.infected = init.individuals,
                             pres.time = 1,
                             table.hosts = iniTable(init.individuals, NA, prefix.host, param.pExit, param.pMove = NA, param.timeContact, param.pTrans, param.moveDist = NA),
                             table.state = NA,
                             type = "singleNone")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)
    if (any(active.hosts)) {

      fun <- function(z) {
        pExitParsed$vect(prestime = pres.time, z[, pExitParsed$vectArgs, with = FALSE])
      }
      p.exit.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][["V1"]]

      exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
    }
    # }

    exiting.full <- active.hosts
    exiting.full[exiting.full] <- exiting

    res$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                    active = 0)]

    active.hosts[active.hosts] <- !exiting # Update active hosts

    if (!any(active.hosts)) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------
    df.meetTransmit <- res$table.hosts[active.hosts, c("hosts.ID")]
    df.meetTransmit[, active.hosts:=hosts.ID]

    fun <- function(z) {
      timeContactParsed$vect(prestime = pres.time, z[, timeContactParsed$vectArgs, with = FALSE])
    }
    timeContact.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

    df.meetTransmit$number.contacts <- timeContact.values

    haveContact <- df.meetTransmit[["number.contacts"]] > 0
    df.meetTransmit <- df.meetTransmit[haveContact]
    active.hosts[active.hosts] <- haveContact # Update active hosts

    if (nrow(df.meetTransmit) > 0) {

      fun <- function(z) {
        pTransParsed$vect(prestime = pres.time, z[, pTransParsed$vectArgs, with = FALSE])
      }

      df.meetTransmit[, "Ptransmit"] <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"] #adds transmission probability to events
      df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Ptransmit"]] > 0] #discards event with probability 0

      if (nrow(df.meetTransmit) > 0) {

        df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

        df.meetTransmit[,"Trans"] <- drawBernouilli(df.meetTransmit[["Ptransmit"]]) #Draws K bernouillis with various probability (see function for more detail)

        df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Trans"]]] #Discards events with no realisation

        if (nrow(df.meetTransmit) >0) {
          table.temp <- vector("list", nrow(df.meetTransmit))
          for (i in 1:nrow(df.meetTransmit)) {

            res$N.infected <- res$N.infected + 1
            hosts.ID <- as.character(paste(prefix.host,res$N.infected,sep="-"))

            table.temp[[i]] <- newLine(hosts.ID, as.character(df.meetTransmit[i,]$active.hosts), NA, pres.time, param.pExit, param.pMove=NA,param.timeContact, param.pTrans,param.moveDist=NA)
          }

          res$table.hosts <- data.table::rbindlist(c(list(res$table.hosts),table.temp))
          data.table::setkey(res$table.hosts,hosts.ID)
        }
      }
    }

    if (progress.bar == TRUE) progressMessage(res$N.infected, pres.time, print.step, length.sim, max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(res$N.infected, pres.time)

  res$pres.time <- pres.time

  return(res)
}
