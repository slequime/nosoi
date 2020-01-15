#' @title Dual-host pathogen in homogeneous hosts populations
#'
#' @description This function, that can be wrapped within \code{\link{nosoiSim}}, runs a dual-host transmission chain simulation, without any structure features in both hosts populations. The simulation stops either at
#' the end of given time (specified by \code{length.sim}) or when the number of hosts infected threshold (\code{max.infected}) is crossed.
#'
#' @section Suffixes:
#' The suffix \code{.A} or \code{.B} specifies if the considered function or parameter concerns host type A or B.
#'
#' @inheritParams singleNone
#' @param max.infected.A specifies the maximum number of individual hosts A that can be infected in the simulation.
#' @param max.infected.B specifies the maximum number of individual hosts B that can be infected in the simulation.
#' @param init.individuals.A number of initially infected individuals (hosts A).
#' @param init.individuals.B number of initially infected individuals (hosts B).
#' @param pExit.A function that gives the probability to exit the simulation for an infected host A (either moving out, dying, etc.).
#' @param param.pExit.A parameter names (list of functions) for the pExit for host-type A.
#' @param timeDep.pExit.A is pExit of host-type A dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param nContact.A function that gives the number of potential transmission events per unit of time for host-type A.
#' @param param.nContact.A parameter names (list of functions) for param.nContact for host-type A.
#' @param timeDep.nContact.A is nContact of host-type A dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param pTrans.A function that gives the probability of transmit a pathogen as a function of time since infection for host A.
#' @param param.pTrans.A parameter names (list of functions) for the pExit  for host A.
#' @param timeDep.pTrans.A is pTrans of host-type A dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param prefix.host.A character(s) to be used as a prefix for the host A identification number.
#' @param pExit.B function that gives the probability to exit the simulation for an infected host B (either moving out, dying, etc.).
#' @param param.pExit.B parameter names (list of functions) for the pExit for host-type B.
#' @param timeDep.pExit.B is pExit of host-type B dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param nContact.B function that gives the number of potential transmission events per unit of time for host B.
#' @param param.nContact.B parameter names (list of functions) for param.nContact for host-type B.
#' @param timeDep.nContact.B is nContact of host-type B dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param pTrans.B function that gives the probability of transmit a pathogen as a function of time since infection for host B.
#' @param param.pTrans.B parameter names (list of functions) for the pExit for host-type B.
#' @param timeDep.pTrans.B is pTrans of host-type B dependent on the absolute time of the simulation (TRUE/FALSE)?
#' @param prefix.host.B character(s) to be used as a prefix for the host B identification number.
#'
#' @inherit singleNone return details
#' @inheritSection singleNone Order of Arguments
#'
#' @seealso For simulations with a discrete structured host population, see \code{\link{dualDiscrete}}. For simulations with a structured population in continuous space, \code{\link{dualContinuous}}
#'
#' @examples
#' \donttest{
#'  #Host A
#' t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
#' pTrans_hostA <- function(t,t_infectA){
#'   if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
#'   if(t/t_infectA > 1){p=0}
#'   return(p)
#' }
#'
#' p_Exit_fctA  <- function(t,t_infectA){
#'   if(t/t_infectA <= 1){p=0}
#'   if(t/t_infectA > 1){p=1}
#'   return(p)
#' }
#'
#' time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}
#'
#' #Host B
#' t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
#' p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#'
#' p_Exit_fct_B  <- function(t,prestime){(sin(prestime/12)+1)/5}
#'
#' pTrans_hostB <- function(t,p_max,t_incub){
#'   if(t <= t_incub){p=0}
#'   if(t >= t_incub){p=p_max}
#'   return(p)
#' }
#'
#' time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}
#'
#' set.seed(90)
#' test.nosoi <- nosoiSim(type="dual", popStructure="none",
#'                        length.sim=40,
#'                        max.infected.A=100,
#'                        max.infected.B=200,
#'                        init.individuals.A=1,
#'                        init.individuals.B=0,
#'                        pExit.A = p_Exit_fctA,
#'                        param.pExit.A = list(t_infectA = t_infectA_fct),
#'                        timeDep.pExit.A=FALSE,
#'                        nContact.A = time_contact_A,
#'                        param.nContact.A = NA,
#'                        timeDep.nContact.A=FALSE,
#'                        pTrans.A = pTrans_hostA,
#'                        param.pTrans.A = list(t_infectA=t_infectA_fct),
#'                                              timeDep.pTrans.A=FALSE,
#'                        prefix.host.A="H",
#'                        pExit.B = p_Exit_fct_B,
#'                        param.pExit.B = NA,
#'                        timeDep.pExit.B=TRUE,
#'                        nContact.B = time_contact_B,
#'                        param.nContact.B = NA,
#'                        timeDep.nContact.B=FALSE,
#'                        pTrans.B = pTrans_hostB,
#'                        param.pTrans.B = list(p_max=p_max_fct_B,
#'                                             t_incub=t_incub_fct_B),
#'                        timeDep.pTrans.B=FALSE,
#'                        prefix.host.B="V")
#'
#'test.nosoi
#' }
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

                     print.progress=TRUE,
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
                             type = "dual",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals.A,
                               table.hosts = iniTable(init.individuals.A, NA, prefix.host.A, ParamHost.A),
                               table.state = NA,
                               prefix.host = prefix.host.A,
                               popStructure = "none"),
                             pop.B = nosoiSimOneConstructor(
                               N.infected = init.individuals.B,
                               table.hosts = iniTable(init.individuals.B, NA, prefix.host.B, ParamHost.B),
                               table.state = NA,
                               prefix.host = prefix.host.B,
                               popStructure = "none"))

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

    if (!any(c(res$host.info.A$table.hosts[["active"]], res$host.info.B$table.hosts[["active"]]))) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------

    #Transmission from A to B
    df.meetTransmit.A <- meetTransmit(res$host.info.A, pres.time, positions = NULL, nContactParsed.A, pTransParsed.A)
    res$host.info.B <- writeInfected(df.meetTransmit.A, res$host.info.B, pres.time, ParamHost.B)

    #Transmission from B to A
    df.meetTransmit.B <- meetTransmit(res$host.info.B, pres.time, positions = NULL, nContactParsed.B, pTransParsed.B)
    res$host.info.A <- writeInfected(df.meetTransmit.B, res$host.info.A, pres.time, ParamHost.A)

    if (print.progress == TRUE) progressMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time=pres.time, print.step=print.step, length.sim=length.sim, max.infected.A=max.infected.A, max.infected.B=max.infected.B, type="dual")
    if (res$host.info.A$N.infected > max.infected.A || res$host.info.B$N.infected > max.infected.B) {break}
  }

  endMessage(Host.count.A=res$host.info.A$N.infected, Host.count.B=res$host.info.B$N.infected, pres.time, type="dual")

  res$total.time <- pres.time

  return(res)
}
