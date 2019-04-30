#' Dual-host without structured host population
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
#' @export dualNone

dualNone <- function(length.sim,
                     max.infected.A,
                     max.infected.B,
                     init.individuals.A,
                     init.individuals.B,

                     pExit.A,
                     param.pExit.A,
                     timeDep.pExit.A=FALSE,
                     timeContact.A,
                     param.timeContact.A,
                     timeDep.timeContact.A=FALSE,
                     pTrans.A,
                     param.pTrans.A,
                     timeDep.pTrans.A=FALSE,
                     prefix.host.A="H",

                     pExit.B,
                     param.pExit.B,
                     timeDep.pExit.B=FALSE,
                     timeContact.B,
                     param.timeContact.B,
                     timeDep.timeContact.B=FALSE,
                     pTrans.B,
                     param.pTrans.B,
                     timeDep.pTrans.B=FALSE,
                     prefix.host.B="V",

                     progress.bar=TRUE,
                     print.step=10,
                     ...){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

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
                             total.time = 1,
                             table.hosts = iniTable(init.individuals, NA, prefix.host, param.pExit, param.pMove = NA, param.timeContact, param.pTrans, param.moveDist = NA),
                             table.state = NA,
                             type = "singleNone")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExitingMoving(res, pres.time, pExitParsed)

    res$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                    active = 0)]

    if (all(res$table.hosts[["active"]] == 0)) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------

    res <- meetTransmit(res,
                        pres.time,
                        positions = NULL,
                        timeContactParsed, pTransParsed,
                        prefix.host, param.pExit, param.pMove = NA, param.timeContact, param.pTrans,
                        param.moveDist = NA)

    if (progress.bar == TRUE) progressMessage(res$N.infected, pres.time, print.step, length.sim, max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(res$N.infected, pres.time)

  res$total.time <- pres.time

  return(res)
}
