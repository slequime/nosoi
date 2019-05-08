#' Dual-host with structured (discrete) host population
#'
#' @description This function runs a dual-host transmission chain simulation, with discrete structure features. The simulation stops either at
#' the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected) is passed.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected.A specifies the maximum number of hosts A that can be infected in the simulation.
#' @param max.infected.B specifies the maximum number of hosts B that can be infected in the simulation.
#' @param init.individuals.A number of initially infected individuals (hosts A).
#' @param init.individuals.B number of initially infected individuals (hosts B).
#' @param init.structure.A which state (e.g. location) the initially infected individuals of host A are located. (NA if init.individual.A is 0)
#' @param init.structure.B which state (e.g. location) the initially infected individuals of host B are located. (NA if init.individual.B is 0)
#' @param structure.raster.A raster object defining the environmental variable for host A.
#' @param structure.raster.B raster object defining the environmental variable for host B.

#' @param pExit.A function that gives the probability to exit the simulation for an infected host A (either moving out, dying, etc.).
#' @param param.pExit.A parameter names (list of functions) for the pExit for host A.
#' @param timeDep.pExit.A is pExit dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param diff.pExit.A is pExit different between states of the structured population (TRUE/FALSE) for host A.
#' @param pMove.A function that gives the probability of a host moving as a function of time for host A.
#' @param param.pMove.A parameter names (list of functions) for the pMove for host A.
#' @param timeDep.pMove.A is pMove dependant on the absolute time of the simulation (TRUE/FALSE) for host A.
#' @param diff.pMove.A is pMove different between states of the structured population (TRUE/FALSE) for host A.
#' @param sdMove.A function that gives the distance travelled (based on coordinates); gives the sd value for the brownian motion for host A.
#' @param param.sdMove.A parameter names (list of functions) for sdMove for host A
#' @param diff.sdMove.A is sdMove dependant on the environmental value (TRUE/FALSE) for host A.
#' @param timeDep.sdMove.A is sdMove dependant on the absolute time of the simulation (TRUE/FALSE) for host A.
#' @param attracted.by.raster.A should the hosts A be attracted by high values in the environmental raster? (TRUE/FALSE)
#' @param nContact.A function that gives the number of potential transmission events per unit of time  for host A.
#' @param param.nContact.A parameter names (list of functions) for param.nContact  for host A.
#' @param timeDep.nContact.A is nContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param diff.nContact.A is nContact different between states of the structured population (TRUE/FALSE) for host A.
#' @param pTrans.A function that gives the probability of transmit a pathogen as a function of time since infection  for host A.
#' @param param.pTrans.A parameter names (list of functions) for the pExit  for host A.
#' @param timeDep.pTrans.A is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)  for host A.
#' @param diff.pTrans.A is pTrans different between states of the structured population (TRUE/FALSE) for host A.
#' @param prefix.host.A character(s) to be used as a prefix for the host A identification number.
#'
#' @param pExit.B function that gives the probability to exit the simulation for an infected host B (either moving out, dying, etc.).
#' @param param.pExit.B parameter names (list of functions) for the pExit for host B.
#' @param timeDep.pExit.B is pExit dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param diff.pExit.B is pExit different between states of the structured population (TRUE/FALSE) for host B.
#' @param pMove.B function that gives the probability of a host moving as a function of time for host B.
#' @param param.pMove.B parameter names (list of functions) for the pMove for host B.
#' @param timeDep.pMove.B is pMove dependant on the absolute time of the simulation (TRUE/FALSE) for host B.
#' @param sdMove.B function that gives the distance travelled (based on coordinates); gives the sd value for the brownian motion for host A.
#' @param param.sdMove.B parameter names (list of functions) for sdMove for host A
#' @param diff.sdMove.B is sdMove dependant on the environmental value (TRUE/FALSE) for host A.
#' @param timeDep.sdMove.B is sdMove dependant on the absolute time of the simulation (TRUE/FALSE) for host A.
#' @param attracted.by.raster.B should the hosts A be attracted by high values in the environmental raster? (TRUE/FALSE)
#' @param diff.pMove.B is pMove different between states of the structured population (TRUE/FALSE) for host B.
#' @param nContact.B function that gives the number of potential transmission events per unit of time  for host B.
#' @param param.nContact.B parameter names (list of functions) for param.nContact  for host B.
#' @param timeDep.nContact.B is nContact dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param diff.nContact.B is nContact different between states of the structured population (TRUE/FALSE) for host B.
#' @param pTrans.B function that gives the probability of transmit a pathogen as a function of time since infection  for host B.
#' @param param.pTrans.B parameter names (list of functions) for the pExit  for host B.
#' @param timeDep.pTrans.B is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)  for host B.
#' @param diff.pTrans.B is pTrans different between states of the structured population (TRUE/FALSE) for host B.
#' @param prefix.host.B character(s) to be used as a prefix for the host B identification number.
#'
#' @param progress.bar if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
#'
#' @export dualContinuous

dualContinuous <- function(length.sim,
                           max.infected.A,
                           max.infected.B,
                           init.individuals.A,
                           init.individuals.B,
                           init.structure.A,
                           init.structure.B,
                           structure.raster.A,
                           structure.raster.B,

                           pExit.A,
                           param.pExit.A,
                           timeDep.pExit.A=FALSE,
                           diff.pExit.A=FALSE,
                           pMove.A,
                           param.pMove.A,
                           timeDep.pMove.A=FALSE,
                           diff.pMove.A=FALSE,
                           sdMove.A,
                           param.sdMove.A,
                           diff.sdMove.A=FALSE,
                           timeDep.sdMove.A=FALSE,
                           attracted.by.raster.A=FALSE,
                           nContact.A,
                           param.nContact.A,
                           timeDep.nContact.A=FALSE,
                           diff.nContact.A=FALSE,
                           pTrans.A,
                           param.pTrans.A,
                           timeDep.pTrans.A=FALSE,
                           diff.pTrans.A=FALSE,
                           prefix.host.A="H",

                           pExit.B,
                           param.pExit.B,
                           timeDep.pExit.B=FALSE,
                           diff.pExit.B=FALSE,
                           pMove.B,
                           param.pMove.B,
                           timeDep.pMove.B=FALSE,
                           diff.pMove.B=FALSE,
                           sdMove.B,
                           param.sdMove.B,
                           diff.sdMove.B=FALSE,
                           timeDep.sdMove.B=FALSE,
                           attracted.by.raster.B=FALSE,
                           nContact.B,
                           param.nContact.B,
                           timeDep.nContact.B=FALSE,
                           diff.nContact.B=FALSE,
                           pTrans.B,
                           param.pTrans.B,
                           timeDep.pTrans.B=FALSE,
                           diff.pTrans.B=FALSE,
                           prefix.host.B="V",

                           progress.bar=TRUE,
                           print.step=10){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected=(max.infected.A+max.infected.B), init.individuals=(init.individuals.A+init.individuals.B))
  none.at.start.A = (init.individuals.A == 0)
  none.at.start.B = (init.individuals.B == 0)

  if(none.at.start.A) init.structure.A = c(0,0)
  if(none.at.start.B) init.structure.B = c(0,0)

  if((!is.function(pMove.A)) && (!is.function(pMove.B))) stop("At least one host must move.")

  #Parsing nContact
  nContactParsed.A <- parseFunction(nContact.A, param.nContact.A, as.character(quote(nContact.A)), diff=diff.nContact.A, timeDep = timeDep.nContact.A, continuous=TRUE)
  nContactParsed.B <- parseFunction(nContact.B, param.nContact.B, as.character(quote(nContact.B)), diff=diff.nContact.B, timeDep = timeDep.nContact.B, continuous=TRUE)

  #Parsing pTrans
  pTransParsed.A <- parseFunction(pTrans.A, param.pTrans.A, as.character(quote(pTrans.A)), diff=diff.pTrans.A, timeDep = timeDep.pTrans.A, continuous=TRUE)
  pTransParsed.B <- parseFunction(pTrans.B, param.pTrans.B, as.character(quote(pTrans.B)), diff=diff.pTrans.B, timeDep = timeDep.pTrans.B, continuous=TRUE)

  #Parsing pExit
  pExitParsed.A <- parseFunction(pExit.A, param.pExit.A, as.character(quote(pExit.A)), diff=diff.pExit.A, timeDep = timeDep.pExit.A, continuous=TRUE)
  pExitParsed.B <- parseFunction(pExit.B, param.pExit.B, as.character(quote(pExit.B)), diff=diff.pExit.B, timeDep = timeDep.pExit.B, continuous=TRUE)

  #Continuous move sanity checks -------------------------------------------------------------------------------------------------------------------

  #Extract environmental value at origin:
  RasterSanityChecks(structure.raster.A,init.structure.A, none.at.start.A)
  RasterSanityChecks(structure.raster.B,init.structure.B, none.at.start.B)

  if(!none.at.start.A) start.env.A <- raster::extract(structure.raster.A,cbind(init.structure.A[1],init.structure.A[2]))
  if(none.at.start.A) start.env.A <- NA
  max.raster.A <- max(structure.raster.A[], na.rm=T)

  if(!none.at.start.B) start.env.B <- raster::extract(structure.raster.B,cbind(init.structure.B[1],init.structure.B[2]))
  if(none.at.start.B) start.env.B <- NA
  max.raster.B <- max(structure.raster.B[], na.rm=T)

  #Parse pMove (same as pExit !!attention if diff)
  if(is.function(pMove.A)) pMoveParsed.A <- parseFunction(pMove.A, param.pMove.A, as.character(quote(pMove.A)), diff=diff.pMove.A, timeDep = timeDep.pMove.A, continuous=TRUE)
  if(is.function(pMove.B)) pMoveParsed.B <- parseFunction(pMove.B, param.pMove.B, as.character(quote(pMove.B)), diff=diff.pMove.B, timeDep = timeDep.pMove.B, continuous=TRUE)

  #Parsing sdMove
  if(is.function(pMove.A)) sdMoveParsed.A <- parseFunction(sdMove.A, param.sdMove.A, as.character(quote(sdMove.A)),diff=diff.sdMove.A, timeDep = timeDep.sdMove.A, continuous=TRUE)
  if(is.function(pMove.B)) sdMoveParsed.B <- parseFunction(sdMove.B, param.sdMove.B, as.character(quote(sdMove.B)),diff=diff.sdMove.B, timeDep = timeDep.sdMove.B, continuous=TRUE)

  #Parsing all parameters
  ParamHost.A <- paramConstructor(param.pExit.A, param.pMove=param.pMove.A, param.nContact.A, param.pTrans.A, param.sdMove=param.sdMove.A)
  ParamHost.B <- paramConstructor(param.pExit.B, param.pMove=param.pMove.B, param.nContact.B, param.pTrans.B, param.sdMove=param.sdMove.B)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             popStructure = "dual",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals.A,
                               table.hosts = iniTable(init.individuals.A, init.structure.A, prefix.host.A, ParamHost.A, current.environmental.value = start.env.A),
                               table.state = iniTableState(init.individuals.A, init.structure.A, prefix.host.A, current.environmental.value = start.env.A),
                               prefix.host = prefix.host.A,
                               geoStructure = "continuous"),
                             pop.B = nosoiSimOneConstructor(
                               N.infected = init.individuals.B,
                               table.hosts = iniTable(init.individuals.B, init.structure.B, prefix.host.B, ParamHost.B, current.environmental.value = start.env.B),
                               table.state = iniTableState(init.individuals.B, init.structure.B, prefix.host.B, current.environmental.value = start.env.B),
                               prefix.host = prefix.host.B,
                               geoStructure = "continuous"))

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

    res$host.info.A <- updateTableState(res$host.info.A, exiting.full.A, pres.time)
    res$host.info.B <- updateTableState(res$host.info.B, exiting.full.B, pres.time)

    if (all(c((res$host.info.A$table.hosts[["active"]] == 0),(res$host.info.B$table.hosts[["active"]] == 0)))) {break}

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    if(is.function(pMove.A)) moving.full.A <- getExitingMoving(res$host.info.A, pres.time, pMoveParsed.A)
    if(is.function(pMove.B)) moving.full.B <- getExitingMoving(res$host.info.B, pres.time, pMoveParsed.B)

    #step 1.2 if moving, where are they going?

    if(is.function(pMove.A)) res$host.info.A <- makeMoves(res$host.info.A, pres.time, moving.full.A,
                                                          sdMoveParsed = sdMoveParsed.A,
                                                          structure.raster = structure.raster.A,
                                                          attracted.by.raster = attracted.by.raster.A,
                                                          max.raster = max.raster.A)

    if(is.function(pMove.B)) res$host.info.B <- makeMoves(res$host.info.B, pres.time, moving.full.B,
                                                          sdMoveParsed = sdMoveParsed.B,
                                                          structure.raster = structure.raster.B,
                                                          attracted.by.raster = attracted.by.raster.B,
                                                          max.raster = max.raster.B)

    #Step 2: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res$host.info.A, pres.time, positions = c("current.in.x", "current.in.y", "current.env.value"), nContactParsed.A, pTransParsed.A)
    res$host.info.B <- writeInfected(df.meetTransmit.A, res$host.info.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res$host.info.B, pres.time, positions = c("current.in.x", "current.in.y", "current.env.value"), nContactParsed.B, pTransParsed.B)
    res$host.info.A <- writeInfected(df.meetTransmit.B, res$host.info.A, pres.time, ParamHost.A)

    if (progress.bar == TRUE) progressMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time=pres.time, print.step=print.step, length.sim=length.sim, max.infected.A=max.infected.A, max.infected.B=max.infected.B, type="dual")
    if (res$host.info.A$N.infected > max.infected.A || res$host.info.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time, type="dual")

  res$total.time <- pres.time

  return(res)
}
