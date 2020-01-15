#' @title Dual-host pathogen in structured (discrete) hosts populations
#'
#' @description This function, that can be wrapped within \code{\link{nosoiSim}}, runs a dual-host transmission chain simulation, with discrete hosts populations structures (e.g. spatial, socio-economic, etc.). The simulation stops either at
#' the end of given time (specified by \code{length.sim}) or when the number of hosts infected threshold (\code{max.infected}) is crossed.
#'
#' @section Structure Matrix:
#' The structure/transition matrix provided provided should of class \code{matrix}, with the same number of rows and columns, rows representing departure state and column the arrival state. All rows should add to 1. Probability values can be different for hosts A and B (so two different matrices), but the name of the column and the rows should be shared.
#' @section Order of Arguments:
#' The user specified function's arguments should follow this order: \code{t} (mandatory), \code{prestime} (optional, only if timeDep is TRUE),
#' \code{current.in} (optional, only if diff is TRUE), \code{host.count.A} or \code{host.count.B} (optional, only if hostCount is TRUE) and \code{parameters} specified in the list.
#'
#' @inheritParams dualNone
#' @param init.structure.A in which state (e.g. location) the initially infected individuals of host-type A are located (\code{NA} if init.individual.A is 0)?
#' @param init.structure.B in which state (e.g. location) the initially infected individuals of host-type B are located (\code{NA} if init.individual.B is 0)?
#' @param structure.matrix.A transition matrix (probabilities) to go from location A (row) to B (column) for host-type A.
#' @param structure.matrix.B transition matrix (probabilities) to go from location A (row) to B (column) for host-type B.
#' @param diff.pExit.A is pExit of host-type A different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pExit.A does pExit of host-type A vary with the host count (of either host-type A or B) in the state? (TRUE/FALSE); diff.pExit.A should be TRUE.
#' @param pMove.A function that gives the probability of a host moving as a function of time for host-type A.
#' @param param.pMove.A parameter names (list of functions) for the pMove for host-type A.
#' @param timeDep.pMove.A is pMove of host-type A dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param diff.pMove.A is pMove of host-type A different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pMove.A does pMove of host-type A vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.pMove.A should be TRUE.
#' @param diff.nContact.A is nContact of host-type A different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.nContact.A does nContact of host-type A vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.nContact.A should be TRUE.
#' @param diff.pTrans.A is pTrans of host-type A different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pTrans.A does pTrans of host-type A vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.pTrans.A should be TRUE.
#' @param diff.pExit.B is pExit of host-type B different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pExit.B does pExit of host-type B vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.pExit.B should be TRUE.
#' @param pMove.B function that gives the probability of a host moving as a function of time for host-type B.
#' @param param.pMove.B parameter names (list of functions) for the pMove for host-type B.
#' @param timeDep.pMove.B is sdMove of host-type B dependent on the absolute time of the simulation (TRUE/FALSE) for host-type B.
#' @param diff.pMove.B is pMove of host-type B different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pMove.B does pMove of host-type B vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.pMove.B should be TRUE.
#' @param diff.nContact.B is nContact of host-type B different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.nContact.B does nContact of host-type B vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.nContact.B should be TRUE.
#' @param diff.pTrans.B is pTrans host-type B different between states of the structured population (TRUE/FALSE)?
#' @param hostCount.pTrans.B does pTrans of host-type B vary with the host count (of either host A or B) in the state? (TRUE/FALSE); diff.pTrans.B should be TRUE.
#'
#' @inherit singleNone return details
#' @inheritSection singleDiscrete Structure Parameters
#' @inheritSection dualNone Suffixes
#'
#' @seealso For simulations with a structure in continuous space, see \code{\link{dualContinuous}}. For simulations without any structures, see \code{\link{dualNone}}.
#'
#' @examples
#' \donttest{
#'#Host A
#'t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
#'pTrans_hostA <- function(t,t_infectA){
#'  if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
#'  if(t/t_infectA > 1){p=0}
#'  return(p)
#'}
#'
#'p_Move_fctA  <- function(t){return(0.1)}
#'
#'p_Exit_fctA  <- function(t,t_infectA){
#'  if(t/t_infectA <= 1){p=0}
#'  if(t/t_infectA > 1){p=1}
#'  return(p)
#'}
#'
#'time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}
#'
#Host B
#'t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
#'p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#'
#'p_Exit_fct_B  <- function(t,current.in){
#'  if(current.in=="A"){return(0.1)}
#'  if(current.in=="B"){return(0.2)}
#'  if(current.in=="C"){return(1)}}
#'
#'pTrans_hostB <- function(t,p_max,t_incub){
#'  if(t <= t_incub){p=0}
#'  if(t >= t_incub){p=p_max}
#'  return(p)
#'}
#'
#'time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}
#'
#'transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),
#'                           nrow = 3, ncol = 3,
#'                           dimnames=list(c("A","B","C"),c("A","B","C")))
#'
#'set.seed(6262)
#'test.nosoi <- nosoiSim(type="dual", popStructure="discrete",
#'                       length.sim=40,
#'                       max.infected.A=100,
#'                       max.infected.B=200,
#'                       init.individuals.A=1,
#'                       init.individuals.B=0,
#'                       init.structure.A="A",
#'                       init.structure.B=NA,
#'                       structure.matrix.A=transition.matrix,
#'                       structure.matrix.B=transition.matrix,
#'                       pExit.A = p_Exit_fctA,
#'                       param.pExit.A = list(t_infectA = t_infectA_fct),
#'                       pMove.A=p_Move_fctA,
#'                       param.pMove.A=NA,
#'                       timeDep.pMove.A=FALSE,
#'                       diff.pMove.A=FALSE,
#'                       timeDep.pExit.A=FALSE,
#'                       nContact.A = time_contact_A,
#'                       param.nContact.A = NA,
#'                       timeDep.nContact.A=FALSE,
#'                       pTrans.A = pTrans_hostA,
#'                       param.pTrans.A = list(t_infectA=t_infectA_fct),
#'                       timeDep.pTrans.A=FALSE,
#'                       prefix.host.A="H",
#'                       pExit.B = p_Exit_fct_B,
#'                       param.pExit.B = NA,
#'                       timeDep.pExit.B=FALSE,
#'                       diff.pExit.B=TRUE,
#'                       pMove.B=NA,
#'                       param.pMove.B=NA,
#'                       timeDep.pMove.B=FALSE,
#'                       diff.pMove.B=FALSE,
#'                       nContact.B = time_contact_B,
#'                       param.nContact.B = NA,
#'                       timeDep.nContact.B=FALSE,
#'                       pTrans.B = pTrans_hostB,
#'                       param.pTrans.B = list(p_max=p_max_fct_B,
#'                                             t_incub=t_incub_fct_B),
#'                       timeDep.pTrans.B=FALSE,
#'                       prefix.host.B="V")
#'
#'test.nosoi
#' }
#' @export dualDiscrete

dualDiscrete <- function(length.sim,
                         max.infected.A,
                         max.infected.B,
                         init.individuals.A,
                         init.individuals.B,
                         init.structure.A,
                         init.structure.B,
                         structure.matrix.A,
                         structure.matrix.B,

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

  if(none.at.start.A) init.structure.A = "NA"
  if(none.at.start.B) init.structure.B = "NA"

  if((!is.function(pMove.A)) && (!is.function(pMove.B))) stop("At least one host must move.")

  #Parsing nContact
  nContactParsed.A <- parseFunction(nContact.A, param.nContact.A, as.character(quote(nContact.A)), diff=diff.nContact.A, timeDep = timeDep.nContact.A, hostCount=hostCount.nContact.A, stateNames=colnames(structure.matrix.A))
  nContactParsed.B <- parseFunction(nContact.B, param.nContact.B, as.character(quote(nContact.B)), diff=diff.nContact.B, timeDep = timeDep.nContact.B, hostCount=hostCount.nContact.B, stateNames=colnames(structure.matrix.B))

  #Parsing pTrans
  pTransParsed.A <- parseFunction(pTrans.A, param.pTrans.A, as.character(quote(pTrans.A)), diff=diff.pTrans.A, timeDep = timeDep.pTrans.A, hostCount=hostCount.pTrans.A, stateNames=colnames(structure.matrix.A))
  pTransParsed.B <- parseFunction(pTrans.B, param.pTrans.B, as.character(quote(pTrans.B)), diff=diff.pTrans.B, timeDep = timeDep.pTrans.B, hostCount=hostCount.pTrans.B, stateNames=colnames(structure.matrix.B))

  #Parsing pExit
  pExitParsed.A <- parseFunction(pExit.A, param.pExit.A, as.character(quote(pExit.A)), diff=diff.pExit.A, timeDep = timeDep.pExit.A, hostCount=hostCount.pExit.A, stateNames=colnames(structure.matrix.A))
  pExitParsed.B <- parseFunction(pExit.B, param.pExit.B, as.character(quote(pExit.B)), diff=diff.pExit.B, timeDep = timeDep.pExit.B, hostCount=hostCount.pExit.B, stateNames=colnames(structure.matrix.B))

  MatrixSanityChecks(structure.matrix.A,init.structure.A, none.at.start.A)
  MatrixSanityChecks(structure.matrix.B,init.structure.B, none.at.start.B)

  #Parse pMove (same as pExit !!attention if diff)
  if(is.function(pMove.A)) pMoveParsed.A <- parseFunction(pMove.A, param.pMove.A, as.character(quote(pMove.A)), diff=diff.pMove.A, timeDep = timeDep.pMove.A, hostCount=hostCount.pMove.A, stateNames=colnames(structure.matrix.A))
  if(is.function(pMove.B)) pMoveParsed.B <- parseFunction(pMove.B, param.pMove.B, as.character(quote(pMove.B)), diff=diff.pMove.B, timeDep = timeDep.pMove.B, hostCount=hostCount.pMove.B, stateNames=colnames(structure.matrix.B))

  #Parsing all parameters
  ParamHost.A <- paramConstructor(param.pExit.A, param.pMove=param.pMove.A, param.nContact.A, param.pTrans.A, param.sdMove=NA)
  ParamHost.B <- paramConstructor(param.pExit.B, param.pMove=param.pMove.B, param.nContact.B, param.pTrans.B, param.sdMove=NA)

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
                                                      current.count.A = init.individuals.A, current.count.B = init.individuals.B),
                               table.state = iniTableState(init.individuals.A, init.structure.A, prefix.host.A),
                               prefix.host = prefix.host.A,
                               popStructure = "discrete"),
                             pop.B = nosoiSimOneConstructor(
                               N.infected = init.individuals.B,
                               table.hosts = iniTable(init.individuals.B, init.structure.B, prefix.host.B, ParamHost.B,
                                                      current.count.A = init.individuals.A, current.count.B = init.individuals.B),
                               table.state = iniTableState(init.individuals.B, init.structure.B, prefix.host.B),
                               prefix.host = prefix.host.B,
                               popStructure = "discrete"))

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
    if(countingHosts) updateHostCount(res$host.info.A, res.B=res$host.info.B, type="discrete")

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    if(is.function(pMove.A)) moving.full.A <- getExitingMoving(res$host.info.A, pres.time, pMoveParsed.A)
    if(is.function(pMove.B)) moving.full.B <- getExitingMoving(res$host.info.B, pres.time, pMoveParsed.B)

    #step 1.2 if moving, where are they going?

    if(is.function(pMove.A)) res$host.info.A <- makeMoves(res$host.info.A, pres.time, moving.full.A, structure.matrix = structure.matrix.A)
    if(is.function(pMove.B)) res$host.info.B <- makeMoves(res$host.info.B, pres.time, moving.full.B, structure.matrix = structure.matrix.B)

    #step 1.3 - count number of hosts per cells
    if(countingHosts) updateHostCount(res$host.info.A, res.B=res$host.info.B, type="discrete")

    #Step 2: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res$host.info.A, pres.time,
                                      positions = c("current.in", "host.count.A", "host.count.B"),
                                      nContactParsed.A, pTransParsed.A)
    res$host.info.B <- writeInfected(df.meetTransmit.A, res$host.info.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res$host.info.B, pres.time,
                                      positions = c("current.in", "host.count.A", "host.count.B"),
                                      nContactParsed.B, pTransParsed.B)
    res$host.info.A <- writeInfected(df.meetTransmit.B, res$host.info.A, pres.time, ParamHost.A)

    #step 1.3 - count number of hosts per cells
    if(countingHosts) updateHostCount(res$host.info.A, res.B = res$host.info.B, type = "discrete")

    if (print.progress == TRUE) progressMessage(Host.count.A = res$host.info.A$N.infected, Host.count.B = res$host.info.B$N.infected,
                                                pres.time = pres.time, print.step = print.step, length.sim = length.sim,
                                                max.infected.A = max.infected.A, max.infected.B = max.infected.B,
                                                type = "dual")
    if (res$host.info.A$N.infected > max.infected.A || res$host.info.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time, type="dual")

  res$total.time <- pres.time

  return(res)
}
