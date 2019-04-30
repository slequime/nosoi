#' Dual-host without structured host population
#'
#' @description This function runs a single-host transmission chain simulation, without any spatial features. The simulation stops either at
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
#' @param timeContact.A function that gives the number of potential transmission events per unit of time  for host A.
#' @param param.timeContact.A parameter names (list of functions) for param.timeContact  for host A.
#' @param timeDep.timeContact.A is timeContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param pTrans.A function that gives the probability of transmit a pathogen as a function of time since infection  for host A.
#' @param param.pTrans.A parameter names (list of functions) for the pExit  for host A.
#' @param timeDep.pTrans.A is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param prefix.host.A character(s) to be used as a prefix for the host A identification number.
#'
#' @param pExit.B function that gives the probability to exit the simulation for an infected host B (either moving out, dying, etc.).
#' @param param.pExit.B parameter names (list of functions) for the pExit for host B.
#' @param timeDep.pExit.B is pExit dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param timeContact.B function that gives the number of potential transmission events per unit of time  for host B.
#' @param param.timeContact.B parameter names (list of functions) for param.timeContact  for host B.
#' @param timeDep.timeContact.B is timeContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
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
                     print.step=10){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected=(max.infected.A+max.infected.B), init.individuals=(init.individuals.A+init.individuals.B))

  #Parsing timeContact
  timeContactParsed.A <- parseFunction(timeContact.A, param.timeContact.A, as.character(quote(timeContact.A)), timeDep = timeDep.timeContact.A)
  timeContactParsed.B <- parseFunction(timeContact.B, param.timeContact.B, as.character(quote(timeContact.B)), timeDep = timeDep.timeContact.B)

  #Parsing pTrans
  pTransParsed.A <- parseFunction(pTrans.A, param.pTrans.A, as.character(quote(pTrans.A)), timeDep = timeDep.pTrans.A)
  pTransParsed.B <- parseFunction(pTrans.B, param.pTrans.B, as.character(quote(pTrans.B)), timeDep = timeDep.pTrans.B)

  #Parsing pExit
  pExitParsed.A <- parseFunction(pExit.A, param.pExit.A, as.character(quote(pExit.A)), timeDep = timeDep.pExit.A)
  pExitParsed.B <- parseFunction(pExit.B, param.pExit.B, as.character(quote(pExit.B)), timeDep = timeDep.pExit.B)

  #Parsing all parameters
  ParamHost.A <- paramConstructor(param.pExit.A, param.pMove=NA, param.timeContact.A, param.pTrans.A, param.moveDist=NA)
  ParamHost.B <- paramConstructor(param.pExit.B, param.pMove=NA, param.timeContact.B, param.pTrans.B, param.moveDist=NA)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res.A <- nosoiSimConstructor(N.infected = init.individuals.A,
                               total.time = 1,
                               table.hosts = iniTable(init.individuals.A, NA, prefix.host.A, ParamHost.A),
                               table.state = NA,
                               prefix.host = prefix.host.A,
                               type = "dualNone")

  res.B <- nosoiSimConstructor(N.infected = init.individuals.B,
                               total.time = 1,
                               table.hosts = iniTable(init.individuals.B, NA, prefix.host.B, ParamHost.B),
                               table.state = NA,
                               prefix.host = prefix.host.B,
                               type = "dualNone")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full.A <- getExitingMoving(res.A, pres.time, pExitParsed.A)
    exiting.full.B <- getExitingMoving(res.B, pres.time, pExitParsed.B)

    res.A$table.hosts[exiting.full.A, `:=` (out.time = as.numeric(pres.time),
                                            active = 0)]
    res.B$table.hosts[exiting.full.B, `:=` (out.time = as.numeric(pres.time),
                                            active = 0)]

    if (all(c((res.A$table.hosts[["active"]] == 0),(res.B$table.hosts[["active"]] == 0)))) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res.A, pres.time, positions = NULL, timeContactParsed.A, pTransParsed.A)
    res.B <- writeInfected(df.meetTransmit.A, res.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res.B, pres.time, positions = NULL, timeContactParsed.B, pTransParsed.B)
    res.A <- writeInfected(df.meetTransmit.B, res.A, pres.time, ParamHost.A)

    if (progress.bar == TRUE) progressMessage(Host.count.A=res.A$N.infected, Host.count.B=res.B$N.infected, pres.time=pres.time, print.step=print.step, length.sim=length.sim, max.infected.A=max.infected.A, max.infected.B=max.infected.B, type="dual")
    if (res.A$N.infected > max.infected.A | res.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res.A$N.infected, Host.count.B=res.B$N.infected, pres.time, type="dual")

  names(res.A) <- paste(names(res.A),"A",sep="_")
  names(res.B) <- paste(names(res.B),"B",sep="_")
  res <- c(res.A,res.B)

  res$total.time <- pres.time
  res$type <- res$type_A

  res[["total.time_A"]] = NULL
  res[["total.time_B"]] = NULL
  res[["type_A"]] = NULL
  res[["type_B"]] = NULL


  return(res)
}
