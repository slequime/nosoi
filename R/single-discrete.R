#' Single-host with structured (discrete) host population
#'
#' @description This function runs a single-host transmission chain simulation, with a structured host population (such as spatial features).
#' The simulation stops either at the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected)
#' is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param init.structure which state (e.g. location) the initially infected individuals are located.
#' @param structure.matrix transition matrix (probabilities) to go from location A (row) to B (column)
#' @param diff.pMove is pMove different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pMove is pMove dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param hostCount.pMove does pMove varies with the host count in the state? (TRUE/FALSE); diff.pMove should be TRUE.
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.nContact is nContact different between states of the structured population (TRUE/FALSE)
#' @param timeDep.nContact is nContact dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param hostCount.nContact does nContact varies with the host count in the state? (TRUE/FALSE); diff.nContact should be TRUE.
#' @param nContact function that gives the number of potential transmission events per unit of time.
#' @param param.nContact parameter names (list of functions) for nContact.
#' @param diff.pTrans is pTrans different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pTrans is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param hostCount.pTrans does pTrans varies with the host count in the state? (TRUE/FALSE); diff.pTrans should be TRUE.
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pTrans.
#' @param diff.pExit is pExit different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pExit is pExit dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param hostCount.pExit  does pExit varies with the host count in the state? (TRUE/FALSE); diff.pExit should be TRUE.
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param print.progress if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step print.progress is TRUE, step with which the progress message will be printed.
#'
#' @export singleDiscrete

singleDiscrete <- function(length.sim,
                           max.infected,
                           init.individuals,
                           init.structure,
                           structure.matrix,
                           diff.pMove=FALSE,
                           timeDep.pMove=FALSE,
                           hostCount.pMove=FALSE,
                           pMove,
                           param.pMove,
                           diff.nContact=FALSE,
                           timeDep.nContact=FALSE,
                           hostCount.nContact=FALSE,
                           nContact,
                           param.nContact,
                           diff.pTrans=FALSE,
                           timeDep.pTrans=FALSE,
                           hostCount.pTrans=FALSE,
                           pTrans,
                           param.pTrans,
                           diff.pExit=FALSE,
                           timeDep.pExit=FALSE,
                           hostCount.pExit=FALSE,
                           pExit,
                           param.pExit,
                           prefix.host="H",
                           print.progress=TRUE,
                           print.step=10){

  #Sanity checks---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  #Parsing nContact
  nContactParsed <- parseFunction(nContact, param.nContact, as.character(quote(nContact)),diff=diff.nContact, timeDep = timeDep.nContact, hostCount=hostCount.nContact, stateNames=colnames(structure.matrix))

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)), diff=diff.pTrans, timeDep = timeDep.pTrans, hostCount=hostCount.pTrans, stateNames=colnames(structure.matrix))

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)), diff=diff.pExit, timeDep = timeDep.pExit, hostCount=hostCount.pExit, stateNames=colnames(structure.matrix))

  #Discrete states sanity checks -------------------------------------------------------------------------------------------------------------------

  MatrixSanityChecks(structure.matrix,init.structure)

  #Parse pMove (same as pExit !!attention if diff)
  pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=diff.pMove, timeDep = timeDep.pMove, hostCount=hostCount.pMove, stateNames=colnames(structure.matrix))

  #Parsing all parameters
  ParamHost <- paramConstructor(param.pExit, param.pMove, param.nContact, param.pTrans, param.sdMove=NA)

  #Are hosts to be counted?
  countingHosts <- any(c(hostCount.pExit, hostCount.pMove, hostCount.nContact, hostCount.pTrans))

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             type = "single",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals,
                               table.hosts = iniTable(init.individuals, init.structure, prefix.host, ParamHost,
                                                      current.count.A = init.individuals),
                               table.state = iniTableState(init.individuals, init.structure, prefix.host),
                               prefix.host = prefix.host,
                               popStructure = "discrete"))

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExitingMoving(res$host.info.A, pres.time, pExitParsed)

    res$host.info.A$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                                    active = 0)]

    res$host.info.A <- updateTableState(res$host.info.A, exiting.full, pres.time)

    if (all(res$host.info.A$table.hosts[["active"]] == 0)) {break}

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, type="discrete")

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    moving.full <- getExitingMoving(res$host.info.A, pres.time, pMoveParsed)

    #step 1.2 if moving, where are they going?

    res$host.info.A <- makeMoves(res$host.info.A, pres.time, moving.full, structure.matrix = structure.matrix)

    #step 1.3 - update host count
    if(countingHosts) updateHostCount(res$host.info.A, type="discrete")

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    df.meetTransmit <- meetTransmit(res$host.info.A, pres.time, positions = c("current.in", "host.count.A"), nContactParsed, pTransParsed)

    res$host.info.A <- writeInfected(df.meetTransmit, res$host.info.A, pres.time, ParamHost)

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, type="discrete")

    if (print.progress == TRUE) progressMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time, print.step = print.step, length.sim = length.sim, max.infected.A = max.infected)
    if (res$host.info.A$N.infected > max.infected) {break}
  }

  endMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time)

  res$total.time <- pres.time
  return(res)
}
