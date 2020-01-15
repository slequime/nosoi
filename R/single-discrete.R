#' @title Single-host pathogen in a structured (discrete) host population
#'
#' @description This function, that can be wrapped within \code{\link{nosoiSim}}, runs a single-host transmission chain simulation, with a discrete host population structure (e.g. spatial, socio-economic, etc.). The simulation stops either at
#' the end of given time (specified by \code{length.sim}) or when the number of hosts infected threshold (\code{max.infected}) is crossed.
#'
#' @section Structure Matrix:
#' The structure matrix provided provided should of class \code{matrix}, with the same number of rows and columns, rows representing departure state and column the arrival state. All rows should add to 1.
#' @section Structure Parameters:
#' The \code{pMove} function should return a single probability (a number between 0 and 1).
#' @section Structure Parameters:
#' The use of \code{diff} (switch to \code{TRUE}) makes the corresponding function use the argument \code{current.in} (for "currently in"). Your function should in that case give a result for every possible discrete state.
#' @section Structure Parameters:
#' The use of \code{hostCount} (switch to \code{TRUE}) makes the corresponding function use the argument \code{host.count}.
#' @section Order of Arguments:
#' The user specified function's arguments should follow this order: \code{t} (mandatory), \code{prestime} (optional, only if timeDep is TRUE),
#' \code{current.in} (optional, only if diff is TRUE), \code{host.count} (optional, only if hostCount is TRUE) and \code{parameters} specified in the list.
#'
#' @inheritParams singleNone
#' @param init.structure in which state (e.g. location) the initially infected individuals are located.
#' @param structure.matrix transition matrix (probabilities) to go from location A (row) to B (column)
#' @param diff.pMove is pMove different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pMove is pMove dependent on the absolute time of the simulation (TRUE/FALSE)
#' @param hostCount.pMove does pMove varies with the host count in the state? (TRUE/FALSE); diff.pMove should be TRUE.
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.nContact is nContact different between states of the structured population (TRUE/FALSE)
#' @param hostCount.nContact does nContact varies with the host count in the state? (TRUE/FALSE); diff.nContact should be TRUE.
#' @param diff.pTrans is pTrans different between states of the structured population (TRUE/FALSE)
#' @param hostCount.pTrans does pTrans varies with the host count in the state? (TRUE/FALSE); diff.pTrans should be TRUE.
#' @param diff.pExit is pExit different between states of the structured population (TRUE/FALSE)
#' @param hostCount.pExit does pExit varies with the host count in the state? (TRUE/FALSE); diff.pExit should be TRUE.
#'
#' @inherit singleNone return details
#'
#' @seealso For simulations with a structure in continuous space, see \code{\link{singleContinuous}}. For simulations without any structures, see \code{\link{singleNone}}.
#'
#' @examples
#' \donttest{
#' t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#' p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#' p_Exit_fct  <- function(t){return(0.08)}
#' p_Move_fct  <- function(t){return(0.1)}
#'
#' proba <- function(t,p_max,t_incub){
#'  if(t <= t_incub){p=0}
#'  if(t >= t_incub){p=p_max}
#'  return(p)
#' }
#'
#' time_contact = function(t){round(rnorm(1, 3, 1), 0)}
#'
#' transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),
#'                            nrow = 3, ncol = 3,
#'                            dimnames=list(c("A","B","C"),c("A","B","C")))
#'
#' set.seed(805)
#' test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
#'                        length=20,
#'                        max.infected=100,
#'                        init.individuals=1,
#'                        init.structure="A",
#'                        structure.matrix=transition.matrix,
#'                        pMove=p_Move_fct,
#'                        param.pMove=NA,
#'                        nContact=time_contact,
#'                        param.nContact=NA,
#'                        pTrans = proba,
#'                        param.pTrans = list(p_max=p_max_fct,
#'                                            t_incub=t_incub_fct),
#'                       pExit=p_Exit_fct,
#'                       param.pExit=NA)
#'}
#' @export singleDiscrete

singleDiscrete <- function(length.sim,
                           max.infected,
                           init.individuals,
                           init.structure,
                           structure.matrix,
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

    res$host.info.A$table.hosts[exiting.full, `:=` (out.time = pres.time,
                                                    active = FALSE)]

    res$host.info.A <- updateTableState(res$host.info.A, exiting.full, pres.time)

    if (!any(res$host.info.A$table.hosts[["active"]])) {break}

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
