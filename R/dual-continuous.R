#' @title Dual-host pathogen in structured (continuous) hosts populations
#'
#' @description This function runs a dual-host transmission chain simulation, with structured hosts populations (such as spatial features) in a shared continuous space.
#' The simulation stops either at the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected)
#' is passed. The movement of hosts on the continuous space map is a random walk (Brownian motion) that can be modified towards a biased random walk where hosts tend to be attracted to higher values of the environmental variable defined by the raster.
#'
#' @section Raster:
#' The structure raster(s) provided provided should of class \code{raster}. High values of the environmental variable can attract hosts if \code{attracted.by.raster} is TRUE. Raster have to share the same space (i.e. also the same cell size and ID).
#' @section Order of Arguments:
#' The user specified function's arguments should follow this order: \code{t} (mandatory), \code{prestime} (optional, only if timeDep is TRUE),
#' \code{current.env.value} (optional, only if diff is TRUE), \code{host.count.A} or \code{host.count.B} (optional, only if hostCount is TRUE) and \code{parameters} specified in the list.
#'
#' @inheritParams dualNone
#' @inheritParams dualDiscrete
#' @param init.structure.A in which location the initially infected host-A individuals are located. A vector of coordinates in the same coordinate space as the raster (NA if init.individual.A is 0).
#' @param init.structure.B in which location the initially infected host-B individuals are located. A vector of coordinates in the same coordinate space as the raster (NA if init.individual.B is 0).
#' @param structure.raster.A raster object defining the environmental variable for host-type A.
#' @param structure.raster.B raster object defining the environmental variable for host B.
#' @param diff.pExit.A does pExit of host-type A depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pExit.A does pExit of host-type A vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.pExit.A should be TRUE.
#' @param diff.pMove.A does pMove of host-type A depend on the environmental variable (set by the raster) (TRUE/FALSE).A.
#' @param hostCount.pMove.A does pMove of host-type A vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.pMove.A should be TRUE.
#' @param sdMove.A function that gives the distance traveled for host-type A (based on coordinates); output is the standard deviation value for the Brownian motion.
#' @param param.sdMove.A parameter names (list of functions) for sdMove for host-type A.
#' @param diff.sdMove.A does sdMove of host-type A depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param timeDep.sdMove.A is sdMove of host-type A dependent on the absolute time of the simulation (TRUE/FALSE) ?
#' @param hostCount.sdMove.A does sdMove varies with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); diff.sdMove.A should be TRUE.
#' @param attracted.by.raster.A should the host-type A be attracted by higher values in the environmental raster? (TRUE/FALSE).
#' @param diff.nContact.A does nContact of host-type A depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.nContact.A does nContact vary with the host count (of either host-type A or B) in each raster cell?? (TRUE/FALSE); diff.nContact.A should be TRUE.
#' @param diff.pTrans.A does pTrans of host-type A depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pTrans.A does pTrans vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); diff.pTrans.A should be TRUE.
#' @param diff.pExit.B does pExit of host-type B depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pExit.B does pExit of host-type B vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.pExit.B should be TRUE.
#' @param diff.pMove.B does pMove of host-type B depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pMove.B does pMove of host-type B vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.pMove.B should be TRUE.
#' @param sdMove.B function that gives the distance traveled for host-type B (based on coordinates); output is the standard deviation value for the Brownian motion.
#' @param param.sdMove.B parameter names (list of functions) for sdMove for host-type B.
#' @param timeDep.sdMove.B is sdMove of host-type B dependent on the absolute time of the simulation (TRUE/FALSE) ?
#' @param diff.sdMove.B does sdMove of host-type B depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.sdMove.B does sdMove of host-type B vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.sdMove.B should be TRUE.
#' @param attracted.by.raster.B should the host-type B be attracted by higher values in the environmental raster? (TRUE/FALSE)
#' @param diff.nContact.B does nContact of host-type B depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.nContact.B does nContact of host-type B vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.nContact.B should be TRUE.
#' @param diff.pTrans.B does pTrans of host-type B depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pTrans.B does pTrans of host-type B vary with the host count (of either host-type A or B) in each raster cell? (TRUE/FALSE); if TRUE, diff.pTrans.B should be TRUE.
#'
#' @inherit singleNone return details
#' @inheritSection singleContinuous Structure Parameters
#' @inheritSection dualNone Suffixes
#'
#' @seealso For simulations with a discrete structure, see \code{\link{dualDiscrete}}. For simulations without any structures, see \code{\link{dualNone}}.
#'
#'
#' @examples
#' \donttest{
#' library(raster)
#'
#'#Generating a raster for the movement
#'set.seed(860)
#'
#'test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
#'test.raster[] <- runif(10000, -80, 180)
#'test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
#'
#'t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#'p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#'p_Move_fct  <- function(t){return(0.1)}
#'
#'sdMove_fct = function(t,current.env.value){return(100/(current.env.value+1))}
#'
#'p_Exit_fct  <- function(t){return(0.08)}
#'
#'proba <- function(t,p_max,t_incub){
#'  if(t <= t_incub){p=0}
#'  if(t >= t_incub){p=p_max}
#'  return(p)
#'}
#'
#'time_contact = function(t){round(rnorm(1, 3, 1), 0)}
#'
#'start.pos <- c(0,0)
#'
#'set.seed(805)
#'test.nosoi <- nosoiSim(type="dual", popStructure="continuous",
#'                       length.sim=200,
#'                       max.infected.A=500,
#'                       max.infected.B=500,
#'                       init.individuals.A=1,
#'                       init.individuals.B=0,
#'                       init.structure.A=start.pos,
#'                       init.structure.B=NA,
#'                       structure.raster.A=test.raster,
#'                       structure.raster.B=test.raster,
#'                       pExit.A=p_Exit_fct,
#'                       param.pExit.A=NA,
#'                       timeDep.pExit.A=FALSE,
#'                       diff.pExit.A=FALSE,
#'                       pMove.A=p_Move_fct,
#'                       param.pMove.A=NA,
#'                       timeDep.pMove.A=FALSE,
#'                       diff.pMove.A=FALSE,
#'                       diff.sdMove.A=TRUE,
#'                       sdMove.A=sdMove_fct,
#'                       param.sdMove.A=NA,
#'                       attracted.by.raster.A=TRUE,
#'                       nContact.A=time_contact,
#'                       param.nContact.A=NA,
#'                       timeDep.nContact.A=FALSE,
#'                       diff.nContact.A=FALSE,
#'                       pTrans.A=proba,
#'                       param.pTrans.A=list(p_max=p_max_fct,
#'                                           t_incub=t_incub_fct),
#'                       timeDep.pTrans.A=FALSE,
#'                       diff.pTrans.A=FALSE,
#'                       prefix.host.A="H",
#'                       pExit.B=p_Exit_fct,
#'                       param.pExit.B=NA,
#'                       timeDep.pExit.B=FALSE,
#'                       diff.pExit.B=FALSE,
#'                       pMove.B=p_Move_fct,
#'                       param.pMove.B=NA,
#'                       timeDep.pMove.B=FALSE,
#'                       diff.pMove.B=FALSE,
#'                       diff.sdMove.B=TRUE,
#'                       sdMove.B=sdMove_fct,
#'                       param.sdMove.B=NA,
#'                       attracted.by.raster.B=TRUE,
#'                       nContact.B=time_contact,
#'                       param.nContact.B=NA,
#'                       timeDep.nContact.B=FALSE,
#'                       diff.nContact.B=FALSE,
#'                       pTrans.B=proba,
#'                       param.pTrans.B=list(p_max=p_max_fct,
#'                                           t_incub=t_incub_fct),
#'                       timeDep.pTrans.B=FALSE,
#'                       diff.pTrans.B=FALSE,
#'                       prefix.host.B="V")
#'test.nosoi
#'}
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
                           hostCount.pExit.A=FALSE,
                           pMove.A,
                           param.pMove.A,
                           timeDep.pMove.A=FALSE,
                           diff.pMove.A=FALSE,
                           hostCount.pMove.A=FALSE,
                           sdMove.A,
                           param.sdMove.A,
                           diff.sdMove.A=FALSE,
                           timeDep.sdMove.A=FALSE,
                           hostCount.sdMove.A=FALSE,
                           attracted.by.raster.A=FALSE,
                           nContact.A,
                           param.nContact.A,
                           timeDep.nContact.A=FALSE,
                           diff.nContact.A=FALSE,
                           hostCount.nContact.A=FALSE,
                           pTrans.A,
                           param.pTrans.A,
                           timeDep.pTrans.A=FALSE,
                           diff.pTrans.A=FALSE,
                           hostCount.pTrans.A=FALSE,
                           prefix.host.A="H",

                           pExit.B,
                           param.pExit.B,
                           timeDep.pExit.B=FALSE,
                           diff.pExit.B=FALSE,
                           hostCount.pExit.B=FALSE,
                           pMove.B,
                           param.pMove.B,
                           timeDep.pMove.B=FALSE,
                           diff.pMove.B=FALSE,
                           hostCount.pMove.B=FALSE,
                           sdMove.B,
                           param.sdMove.B,
                           diff.sdMove.B=FALSE,
                           timeDep.sdMove.B=FALSE,
                           hostCount.sdMove.B=FALSE,
                           attracted.by.raster.B=FALSE,
                           nContact.B,
                           param.nContact.B,
                           timeDep.nContact.B=FALSE,
                           diff.nContact.B=FALSE,
                           hostCount.nContact.B=FALSE,
                           pTrans.B,
                           param.pTrans.B,
                           timeDep.pTrans.B=FALSE,
                           diff.pTrans.B=FALSE,
                           hostCount.pTrans.B=FALSE,
                           prefix.host.B="V",

                           print.progress=TRUE,
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
  nContactParsed.A <- parseFunction(nContact.A, param.nContact.A, as.character(quote(nContact.A)), diff=diff.nContact.A, timeDep = timeDep.nContact.A, hostCount=hostCount.nContact.A,continuous=TRUE)
  nContactParsed.B <- parseFunction(nContact.B, param.nContact.B, as.character(quote(nContact.B)), diff=diff.nContact.B, timeDep = timeDep.nContact.B, hostCount=hostCount.nContact.B, continuous=TRUE)

  #Parsing pTrans
  pTransParsed.A <- parseFunction(pTrans.A, param.pTrans.A, as.character(quote(pTrans.A)), diff=diff.pTrans.A, timeDep = timeDep.pTrans.A, hostCount=hostCount.pTrans.A, continuous=TRUE)
  pTransParsed.B <- parseFunction(pTrans.B, param.pTrans.B, as.character(quote(pTrans.B)), diff=diff.pTrans.B, timeDep = timeDep.pTrans.B, hostCount=hostCount.pTrans.B, continuous=TRUE)

  #Parsing pExit
  pExitParsed.A <- parseFunction(pExit.A, param.pExit.A, as.character(quote(pExit.A)), diff=diff.pExit.A, timeDep = timeDep.pExit.A, hostCount=hostCount.pExit.A, continuous=TRUE)
  pExitParsed.B <- parseFunction(pExit.B, param.pExit.B, as.character(quote(pExit.B)), diff=diff.pExit.B, timeDep = timeDep.pExit.B, hostCount=hostCount.pExit.B, continuous=TRUE)

  #Continuous move sanity checks -------------------------------------------------------------------------------------------------------------------

  #Extract environmental value at origin:
  RasterSanityChecks(structure.raster.A,init.structure.A, none.at.start.A)
  RasterSanityChecks(structure.raster.B,init.structure.B, none.at.start.B)


  if(!none.at.start.A) {
    start.cell.A <-  raster::cellFromXY(structure.raster.A,cbind(init.structure.A[1],init.structure.A[2]))
    start.env.A <- raster::extract(structure.raster.A,start.cell.A)
  }

  if(none.at.start.A) {
    start.cell.A <- NA
    start.env.A <- NA
  }

  max.raster.A <- max(structure.raster.A[], na.rm=T)

  if(!none.at.start.B) {
    start.cell.B <-  raster::cellFromXY(structure.raster.B,cbind(init.structure.B[1],init.structure.B[2]))
    start.env.B <- raster::extract(structure.raster.B,start.cell.B)
  }

  if(none.at.start.B) {
    start.cell.B <- NA
    start.env.B <- NA
  }

  max.raster.B <- max(structure.raster.B[], na.rm=T)

  #Parse pMove (same as pExit !!attention if diff)
  if(is.function(pMove.A)) pMoveParsed.A <- parseFunction(pMove.A, param.pMove.A, as.character(quote(pMove.A)), diff=diff.pMove.A, timeDep = timeDep.pMove.A, hostCount=hostCount.pMove.A, continuous=TRUE)
  if(is.function(pMove.B)) pMoveParsed.B <- parseFunction(pMove.B, param.pMove.B, as.character(quote(pMove.B)), diff=diff.pMove.B, timeDep = timeDep.pMove.B, hostCount=hostCount.pMove.B, continuous=TRUE)

  #Parsing sdMove
  if(is.function(pMove.A)) sdMoveParsed.A <- parseFunction(sdMove.A, param.sdMove.A, as.character(quote(sdMove.A)),diff=diff.sdMove.A, timeDep = timeDep.sdMove.A, hostCount=hostCount.sdMove.A, continuous=TRUE)
  if(is.function(pMove.B)) sdMoveParsed.B <- parseFunction(sdMove.B, param.sdMove.B, as.character(quote(sdMove.B)),diff=diff.sdMove.B, timeDep = timeDep.sdMove.B, hostCount=hostCount.sdMove.B, continuous=TRUE)

  #Parsing all parameters
  ParamHost.A <- paramConstructor(param.pExit.A, param.pMove=param.pMove.A, param.nContact.A, param.pTrans.A, param.sdMove=param.sdMove.A)
  ParamHost.B <- paramConstructor(param.pExit.B, param.pMove=param.pMove.B, param.nContact.B, param.pTrans.B, param.sdMove=param.sdMove.B)

  #Are hosts to be counted?
  countingHosts <- any(c(hostCount.pExit.A, hostCount.pMove.A, hostCount.nContact.A, hostCount.pTrans.A),c(hostCount.pExit.B, hostCount.pMove.B, hostCount.nContact.B, hostCount.pTrans.B))

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             type = "dual",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals.A,
                               table.hosts = iniTable(init.individuals.A, init.structure.A, prefix.host.A, ParamHost.A,
                                                      current.environmental.value = start.env.A, current.cell.number.raster = start.cell.A,
                                                      current.count.A = init.individuals.A, current.count.B = init.individuals.B),
                               table.state = iniTableState(init.individuals.A, init.structure.A, prefix.host.A,
                                                           current.environmental.value = start.env.A, current.cell.number.raster = start.cell.A),
                               prefix.host = prefix.host.A,
                               popStructure = "continuous"),
                             pop.B = nosoiSimOneConstructor(
                               N.infected = init.individuals.B,
                               table.hosts = iniTable(init.individuals.B, init.structure.B, prefix.host.B, ParamHost.B,
                                                      current.environmental.value = start.env.B, current.cell.number.raster = start.cell.B,
                                                      current.count.A = init.individuals.A, current.count.B = init.individuals.B),
                               table.state = iniTableState(init.individuals.B, init.structure.B, prefix.host.B,
                                                           current.environmental.value = start.env.B, current.cell.number.raster = start.cell.B),
                               prefix.host = prefix.host.B,
                               popStructure = "continuous"))

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full.A <- getExitingMoving(res$host.info.A, pres.time, pExitParsed.A)
    exiting.full.B <- getExitingMoving(res$host.info.B, pres.time, pExitParsed.B)

    res$host.info.A$table.hosts[exiting.full.A, `:=` (out.time = pres.time,
                                                      active = FALSE)]
    res$host.info.B$table.hosts[exiting.full.B, `:=` (out.time = pres.time,
                                                      active = FALSE)]

    res$host.info.A <- updateTableState(res$host.info.A, exiting.full.A, pres.time)
    res$host.info.B <- updateTableState(res$host.info.B, exiting.full.B, pres.time)

    if (!any(c(res$host.info.A$table.hosts[["active"]], res$host.info.B$table.hosts[["active"]]))) {break}

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, res.B=res$host.info.B, type="continuous")

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

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, res.B=res$host.info.B, type="continuous")

    #Step 2: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res$host.info.A, pres.time,
                                      positions = c("current.in.x", "current.in.y",
                                                    "current.env.value","current.cell.raster",
                                                    "host.count.A", "host.count.B"),
                                      nContactParsed.A, pTransParsed.A)
    res$host.info.B <- writeInfected(df.meetTransmit.A, res$host.info.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res$host.info.B, pres.time,
                                      positions = c("current.in.x", "current.in.y",
                                                    "current.env.value","current.cell.raster",
                                                    "host.count.A","host.count.B"),
                                      nContactParsed.B, pTransParsed.B)
    res$host.info.A <- writeInfected(df.meetTransmit.B, res$host.info.A, pres.time, ParamHost.A)

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, res.B=res$host.info.B, type="continuous")

    if (print.progress == TRUE) progressMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time=pres.time, print.step=print.step, length.sim=length.sim, max.infected.A=max.infected.A, max.infected.B=max.infected.B, type="dual")
    if (res$host.info.A$N.infected > max.infected.A || res$host.info.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time, type="dual")

  res$total.time <- pres.time

  return(res)
}
