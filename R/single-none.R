#' Single-host without structured host population
#'
#' @description This function runs a single-host transmission chain simulation, without any spatial features. The simulation stops either at
#' the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected) is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param nContact function that gives the number of potential transmission events per unit of time.
#' @param param.nContact parameter names (list of functions) for param.nContact.
#' @param timeDep.nContact is nContact dependant on the absolute time of the simulation (TRUE/FALSE)
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
                       nContact,
                       param.nContact,
                       timeDep.nContact=FALSE,
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

  #Parsing nContact
  nContactParsed <- parseFunction(nContact, param.nContact, as.character(quote(nContact)),timeDep=timeDep.nContact)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),timeDep=timeDep.pTrans)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),timeDep=timeDep.pExit)

  #Parsing all parameters
  ParamHost <- paramConstructor(param.pExit, param.pMove=NA, param.nContact, param.pTrans, param.coordMove=NA)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(N.infected = init.individuals,
                             total.time = 1,
                             table.hosts = iniTable(init.individuals, NA, prefix.host, ParamHost),
                             table.state = NA,
                             prefix.host = prefix.host,
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

    df.meetTransmit <- meetTransmit(res, pres.time, positions = NULL, nContactParsed, pTransParsed)

    res <- writeInfected(df.meetTransmit, res, pres.time, ParamHost)

    if (progress.bar == TRUE) progressMessage(Host.count.A = res$N.infected, pres.time = pres.time, print.step = print.step, length.sim = length.sim, max.infected.A = max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(Host.count.A = res$N.infected, pres.time = pres.time)

  res$total.time <- pres.time

  return(res)
}
