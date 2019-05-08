#' Dual-host without structured host population
#'
#' @description This function runs a dual-host transmission chain simulation, without any spatial features. The simulation stops either at
#' the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected) is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected.A specifies the maximum number of hosts A that can be infected in the simulation.
#' @param max.infected.B specifies the maximum number of hosts B that can be infected in the simulation.
#' @param init.individuals.A number of initially infected individuals (hosts A).
#' @param init.individuals.B number of initially infected individuals (hosts B).

#' @param pExit.A function that gives the probability to exit the simulation for an infected host A (either moving out, dying, etc.).
#' @param param.pExit.A parameter names (list of functions) for the pExit for host A.
#' @param timeDep.pExit.A is pExit dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param nContact.A function that gives the number of potential transmission events per unit of time  for host A.
#' @param param.nContact.A parameter names (list of functions) for param.nContact  for host A.
#' @param timeDep.nContact.A is nContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param pTrans.A function that gives the probability of transmit a pathogen as a function of time since infection  for host A.
#' @param param.pTrans.A parameter names (list of functions) for the pExit  for host A.
#' @param timeDep.pTrans.A is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param prefix.host.A character(s) to be used as a prefix for the host A identification number.
#'
#' @param pExit.B function that gives the probability to exit the simulation for an infected host B (either moving out, dying, etc.).
#' @param param.pExit.B parameter names (list of functions) for the pExit for host B.
#' @param timeDep.pExit.B is pExit dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param nContact.B function that gives the number of potential transmission events per unit of time  for host B.
#' @param param.nContact.B parameter names (list of functions) for param.nContact  for host B.
#' @param timeDep.nContact.B is nContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param pTrans.B function that gives the probability of transmit a pathogen as a function of time since infection  for host B.
#' @param param.pTrans.B parameter names (list of functions) for the pExit  for host B.
#' @param timeDep.pTrans.B is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param prefix.host.B character(s) to be used as a prefix for the host B identification number.
#'
#' @param progress.bar if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
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
                     nContact.A,
                     param.nContact.A,
                     timeDep.nContact.A=FALSE,
                     pTrans.A,
                     param.pTrans.A,
                     timeDep.pTrans.A=FALSE,
                     prefix.host.A="H",

                     pExit.B,
                     param.pExit.B,
                     timeDep.pExit.B=FALSE,
                     nContact.B,
                     param.nContact.B,
                     timeDep.nContact.B=FALSE,
                     pTrans.B,
                     param.pTrans.B,
                     timeDep.pTrans.B=FALSE,
                     prefix.host.B="V",

                     progress.bar=TRUE,
                     print.step=10){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected=(max.infected.A+max.infected.B), init.individuals=(init.individuals.A+init.individuals.B))

  #Parsing nContact
  nContactParsed.A <- parseFunction(nContact.A, param.nContact.A, as.character(quote(nContact.A)), timeDep = timeDep.nContact.A)
  nContactParsed.B <- parseFunction(nContact.B, param.nContact.B, as.character(quote(nContact.B)), timeDep = timeDep.nContact.B)

  #Parsing pTrans
  pTransParsed.A <- parseFunction(pTrans.A, param.pTrans.A, as.character(quote(pTrans.A)), timeDep = timeDep.pTrans.A)
  pTransParsed.B <- parseFunction(pTrans.B, param.pTrans.B, as.character(quote(pTrans.B)), timeDep = timeDep.pTrans.B)

  #Parsing pExit
  pExitParsed.A <- parseFunction(pExit.A, param.pExit.A, as.character(quote(pExit.A)), timeDep = timeDep.pExit.A)
  pExitParsed.B <- parseFunction(pExit.B, param.pExit.B, as.character(quote(pExit.B)), timeDep = timeDep.pExit.B)

  #Parsing all parameters
  ParamHost.A <- paramConstructor(param.pExit.A, param.pMove=NA, param.nContact.A, param.pTrans.A, param.sdMove=NA)
  ParamHost.B <- paramConstructor(param.pExit.B, param.pMove=NA, param.nContact.B, param.pTrans.B, param.sdMove=NA)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             popStructure = "dual",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals.A,
                               table.hosts = iniTable(init.individuals.A, NA, prefix.host.A, ParamHost.A),
                               table.state = NA,
                               prefix.host = prefix.host.A,
                               geoStructure = "none"),
                             pop.B = nosoiSimOneConstructor(
                               N.infected = init.individuals.B,
                               table.hosts = iniTable(init.individuals.B, NA, prefix.host.B, ParamHost.B),
                               table.state = NA,
                               prefix.host = prefix.host.B,
                               geoStructure = "none"))

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full.A <- getExitingMoving(res$host.info.A, pres.time, pExitParsed.A)
    exiting.full.B <- getExitingMoving(res$host.info.B, pres.time, pExitParsed.B)

    res$host.info.A$table.hosts[exiting.full.A, `:=` (out.time = as.numeric(pres.time),
                                                      active = 0)]
    res$host.info.B$table.hosts[exiting.full.B, `:=` (out.time = as.numeric(pres.time),
                                                      active = 0)]

    if (all(c((res$host.info.A$table.hosts[["active"]] == 0),(res$host.info.B$table.hosts[["active"]] == 0)))) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res$host.info.A, pres.time, positions = NULL, nContactParsed.A, pTransParsed.A)
    res$host.info.B <- writeInfected(df.meetTransmit.A, res$host.info.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res$host.info.B, pres.time, positions = NULL, nContactParsed.B, pTransParsed.B)
    res$host.info.A <- writeInfected(df.meetTransmit.B, res$host.info.A, pres.time, ParamHost.A)

    if (progress.bar == TRUE) progressMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time=pres.time, print.step=print.step, length.sim=length.sim, max.infected.A=max.infected.A, max.infected.B=max.infected.B, type="dual")
    if (res$host.info.A$N.infected > max.infected.A || res$host.info.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time, type="dual")

  res$total.time <- pres.time

  return(res)
}
