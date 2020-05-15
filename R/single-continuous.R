#' @title Single-host pathogen in a structured (continuous) host population
#'
#' @description This function runs a single-host transmission chain simulation, with a structured host population (such as spatial features) in a continuous space.
#' The simulation stops either at the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected)
#' is passed. The movement of hosts on the continuous space map is a random walk (Brownian motion) that can be modified towards a biased random walk where hosts tend to be attracted to higher values of the environmental variable defined by the raster.
#'
#' @section Raster:
#' The structure raster provided provided should of class \code{raster}. High values of the environmental variable can attract hosts if \code{attracted.by.raster} is TRUE.
#' @section Structure Parameters:
#' The \code{pMove} function should return a single probability (a number between 0 and 1), and \code{sdMove} a real number (keep in mind this number is related to your coordinate space).
#' @section Structure Parameters:
#' The use of \code{diff} (switch to \code{TRUE}) makes the corresponding function use the argument \code{current.env.value} (for "current environmental value").
#' @section Structure Parameters:
#' The use of \code{hostCount} (switch to \code{TRUE}) makes the corresponding function use the argument \code{host.count}.
#' @section Order of Arguments:
#' The user specified function's arguments should follow this order: \code{t} (mandatory), \code{prestime} (optional, only if timeDep is TRUE),
#' \code{current.env.value} (optional, only if diff is TRUE), \code{host.count} (optional, only if hostCount is TRUE) and \code{parameters} specified in the list.
#'
#' @inheritParams singleNone
#' @inheritParams singleDiscrete
#' @param init.structure in which location the initially infected individuals are located. A vector of coordinates in the same coordinate space as the raster.
#' @param structure.raster raster object defining the environmental variable.
#' @param diff.pMove does pMove depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param timeDep.pMove does pMove depend on the absolute time of the simulation (TRUE/FALSE).
#' @param hostCount.pMove does pMove vary with the host count in each raster cell? (TRUE/FALSE); if TRUE, diff.pMove should also be TRUE.
#' @param diff.sdMove does sdMove depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param timeDep.sdMove does sdMove depend on the absolute time of the simulation (TRUE/FALSE).
#' @param hostCount.sdMove does sdMove vary with the host count in each raster cell? (TRUE/FALSE); if TRUE, diff.sdMove should be TRUE.
#' @param sdMove function that gives the distance traveled (based on coordinates); output is the standard deviation value for the Brownian motion.
#' @param param.sdMove parameter names (list of functions) for sdMove.
#' @param attracted.by.raster should the hosts be attracted by higher values in the environmental raster? (TRUE/FALSE).
#' @param diff.nContact does nContact depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.nContact does nContact vary with the host count in each raster cell? (TRUE/FALSE); if TRUE, diff.nContact should be TRUE.
#' @param diff.pTrans does pTrans depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pTrans does pTrans vary with the host count in each raster cell? (TRUE/FALSE); if TRUE, diff.pTrans should be TRUE.
#' @param diff.pExit does pExit depend on the environmental variable (set by the raster) (TRUE/FALSE).
#' @param hostCount.pExit does pExit vary with the host count in each raster cell? (TRUE/FALSE); if TRUE, diff.pExit should be TRUE.
#'
#' @inherit singleNone return details
#'
#' @seealso For simulations with a discrete structure, see \code{\link{singleDiscrete}}. For simulations without any structures, see \code{\link{singleNone}}.
#'
#' @examples
#' \donttest{
#' library(raster)
#' #Generating a raster for the movement
#' set.seed(860)
#'
#' test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
#' test.raster[] <- runif(10000, -80, 180)
#' test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
#' plot(test.raster)
#'
#' t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#' p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#' p_Move_fct  <- function(t){return(0.1)}
#'
#' sdMove_fct = function(t,current.env.value){return(100/(current.env.value+1))}
#'
#' p_Exit_fct  <- function(t){return(0.08)}
#'
#' proba <- function(t,p_max,t_incub){
#'   if(t <= t_incub){p=0}
#'   if(t >= t_incub){p=p_max}
#'   return(p)
#' }
#'
#' time_contact = function(t){round(rnorm(1, 3, 1), 0)}
#'
#' start.pos <- c(0,0)
#'
#' test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
#'                    length=200,
#'                    max.infected=500,
#'                    init.individuals=1,
#'                    init.structure=start.pos,
#'                    structure.raster=test.raster,
#'                    pMove=p_Move_fct,
#'                    param.pMove=NA,
#'                    diff.sdMove=TRUE,
#'                    sdMove=sdMove_fct,
#'                    param.sdMove=NA,
#'                    attracted.by.raster=TRUE,
#'                    nContact=time_contact,
#'                    param.nContact=NA,
#'                    pTrans = proba,
#'                    param.pTrans = list(p_max=p_max_fct,
#'                                        t_incub=t_incub_fct),
#'                    pExit=p_Exit_fct,
#'                    param.pExit=NA)
#' }
#' @export singleContinuous

singleContinuous <- function(length.sim,
                             max.infected,
                             init.individuals,
                             init.structure,
                             structure.raster,
                             diff.pExit=FALSE,
                             timeDep.pExit=FALSE,
                             hostCount.pExit=FALSE,
                             pExit,
                             param.pExit,
                             diff.pMove=FALSE,
                             timeDep.pMove=FALSE,
                             hostCount.pMove=FALSE,
                             pMove,
                             param.pMove,
                             diff.sdMove=FALSE,
                             timeDep.sdMove=FALSE,
                             hostCount.sdMove=FALSE,
                             sdMove,
                             param.sdMove,
                             attracted.by.raster=FALSE,
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
                             prefix.host="H",
                             print.progress=TRUE,
                             print.step=10){

  #Sanity checks---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  #Parsing nContact
  nContactParsed <- parseFunction(nContact, param.nContact, as.character(quote(nContact)),diff=diff.nContact, timeDep = timeDep.nContact, hostCount = hostCount.nContact, continuous=TRUE)

  #Parsing sdMove
  sdMoveParsed <- parseFunction(sdMove, param.sdMove, as.character(quote(sdMove)),diff=diff.sdMove, timeDep = timeDep.sdMove, hostCount= hostCount.sdMove, continuous=TRUE)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=diff.pTrans, timeDep = timeDep.pTrans, hostCount= hostCount.pTrans, continuous=TRUE)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=diff.pExit, timeDep = timeDep.pExit, hostCount= hostCount.pExit, continuous=TRUE)

  #Continuous states sanity checks -------------------------------------------------------------------------------------------------------------------

  #Extract environmental value at origin:
  RasterSanityChecks(structure.raster,init.structure)
  start.cell <- raster::cellFromXY(structure.raster,cbind(init.structure[1],init.structure[2]))
  start.env <- raster::extract(structure.raster,start.cell)
  max.raster <- max(structure.raster[], na.rm=T)

  #Parse pMove (same as pExit !!attention if diff)
  pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=diff.pMove, timeDep = timeDep.pMove, hostCount= hostCount.pMove,  continuous=TRUE)

  #Parsing all parameters
  ParamHost <- paramConstructor(param.pExit, param.pMove, param.nContact, param.pTrans, param.sdMove)

  #Are hosts to be counted?
  countingHosts <- any(c(hostCount.pExit, hostCount.pMove, hostCount.sdMove, hostCount.nContact, hostCount.pTrans))

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             type = "single",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals,
                               table.hosts = iniTable(init.individuals, init.structure, prefix.host, ParamHost,
                                                      current.environmental.value = start.env, current.cell.number.raster = start.cell,
                                                      current.count.A = init.individuals),
                               table.state = iniTableState(init.individuals, init.structure, prefix.host,
                                                           current.environmental.value = start.env, current.cell.number.raster = start.cell),
                               prefix.host = prefix.host,
                               popStructure = "continuous"))

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExitingMoving(res$host.info.A, pres.time, pExitParsed)

    res$host.info.A$table.hosts[exiting.full, `:=` (out.time = pres.time,
                                                    active = FALSE)]

    res$host.info.A <- updateTableState(res$host.info.A, exiting.full, pres.time)

    if (!any(res$host.info.A$table.hosts[["active"]])) {break}

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, type="continuous")

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    moving.full <- getExitingMoving(res$host.info.A, pres.time, pMoveParsed)

    #step 1.2 if moving, where are they going?

    res$host.info.A <- makeMoves(res$host.info.A, pres.time, moving.full,
                                 sdMoveParsed = sdMoveParsed,
                                 structure.raster = structure.raster,
                                 attracted.by.raster = attracted.by.raster,
                                 max.raster = max.raster)

    #step 1.3 - count number of hosts per cells
    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, type="continuous")

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    df.meetTransmit <- meetTransmit(res$host.info.A, pres.time,
                                    positions = c("current.in.x", "current.in.y",
                                                  "current.env.value","current.cell.raster",
                                                  "host.count.A"),
                                    nContactParsed, pTransParsed)

    res$host.info.A <- writeInfected(df.meetTransmit, res$host.info.A, pres.time, ParamHost)

    #update host.count
    if(countingHosts) updateHostCount(res$host.info.A, type="continuous")

    if (print.progress == TRUE) progressMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time, print.step = print.step, length.sim = length.sim, max.infected.A = max.infected)
    if (res$host.info.A$N.infected > max.infected) {break}
  }

  endMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time)

  res$total.time <- pres.time

  return(res)
}
