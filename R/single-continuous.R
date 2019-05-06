#' Single-host with structured host population
#'
#' @description This function runs a single-host transmission chain simulation, with a structured host population (such as spatial features).
#' The simulation stops either at the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected)
#' is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param init.structure which state (i.e. location) the initially infected individuals are located.
#' @param structure.raster raster object defining the environmental variable.
#' @param diff.pMove is pMove different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pMove is pMove dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.coordMove is coordMove dependant on the environmental value (TRUE/FALSE).
#' @param timeDep.coordMove is coordMove dependant on the absolute time of the simulation (TRUE/FALSE).
#' @param coordMove function that gives the distance travelled (based on coordinates); gives the sd value for the brownian motion.
#' @param param.coordMove parameter names (list of functions) for coordMove.
#' @param attracted.by.raster should the hosts be attracted by high values in the environmental raster? (TRUE/FALSE)
#' @param diff.nContact is nContact different between states of the structured population (TRUE/FALSE)
#' @param timeDep.nContact is nContact dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param nContact function that gives the number of potential transmission events per unit of time.
#' @param param.nContact parameter names (list of functions) for nContact.
#' @param diff.pTrans is pTrans different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pTrans is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pTrans.
#' @param diff.pExit is pExit different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pExit is pExit dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param progress.bar if TRUE, print message on simulation state to screen.
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @export singleContinuous

singleContinuous <- function(length.sim,
                             max.infected,
                             init.individuals,
                             init.structure,
                             structure.raster,
                             diff.pMove=FALSE,
                             timeDep.pMove=FALSE,
                             pMove,
                             param.pMove,
                             diff.coordMove=FALSE,
                             timeDep.coordMove=FALSE,
                             coordMove,
                             param.coordMove,
                             attracted.by.raster=FALSE,
                             diff.nContact=FALSE,
                             timeDep.nContact=FALSE,
                             nContact,
                             param.nContact,
                             diff.pTrans=FALSE,
                             timeDep.pTrans=FALSE,
                             pTrans,
                             param.pTrans,
                             diff.pExit=FALSE,
                             timeDep.pExit=FALSE,
                             pExit,
                             param.pExit,
                             prefix.host="H",
                             progress.bar=TRUE,
                             print.step=10,
                             ...){

  #Sanity checks---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  #Parsing nContact
  nContactParsed <- parseFunction(nContact, param.nContact, as.character(quote(nContact)),diff=diff.nContact, timeDep = timeDep.nContact, continuous=TRUE)

  #Parsing coordMove
  coordMoveParsed <- parseFunction(coordMove, param.coordMove, as.character(quote(coordMove)),diff=diff.coordMove, timeDep = timeDep.coordMove, continuous=TRUE)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=diff.pTrans, timeDep = timeDep.pTrans, continuous=TRUE)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=diff.pExit, timeDep = timeDep.pExit, continuous=TRUE)

  #Discrete states sanity checks -------------------------------------------------------------------------------------------------------------------

  #Extract environmental value at origin:
  RasterSanityChecks(structure.raster,init.structure)
  start.env <- raster::extract(structure.raster,cbind(init.structure[1],init.structure[2]))
  max.raster <- max(structure.raster[], na.rm=T)

  #Parse pMove (same as pExit !!attention if diff)
  pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=diff.pMove, timeDep = timeDep.pMove, continuous=TRUE)

  #Parsing all parameters
  ParamHost <- paramConstructor(param.pExit, param.pMove, param.nContact, param.pTrans, param.coordMove)

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(N.infected = init.individuals,
                             total.time = 1,
                             table.hosts = iniTable(init.individuals, init.structure, prefix.host, ParamHost, current.environmental.value = start.env),
                             table.state = iniTableState(init.individuals, init.structure, prefix.host, current.environmental.value = start.env),
                             prefix.host = prefix.host,
                             type = "singleContinuous")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExitingMoving(res, pres.time, pExitParsed)

    res$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                        active = 0)]

    res <- updateTableState(res, exiting.full, pres.time)

    if (all(res$table.hosts[["active"]] == 0)) {break}

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    moving.full <- getExitingMoving(res, pres.time, pMoveParsed)

    #step 1.2 if moving, where are they going?

    res <- makeMoves(res, pres.time, moving.full,
                     coordMoveParsed = coordMoveParsed,
                     structure.raster = structure.raster,
                     attracted.by.raster = attracted.by.raster,
                     max.raster = max.raster)

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    df.meetTransmit <- meetTransmit(res, pres.time, positions = c("current.in.x", "current.in.y", "current.env.value"), nContactParsed, pTransParsed)

    res <- writeInfected(df.meetTransmit, res, pres.time, ParamHost)


    if (progress.bar == TRUE) progressMessage(Host.count.A = res$N.infected, pres.time = pres.time, print.step = print.step, length.sim = length.sim, max.infected.A = max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(Host.count.A = res$N.infected, pres.time = pres.time)

  res$total.time <- pres.time

  return(res)
}
