#' @title Single-host pathogen in a homogeneous host population
#'
#' @description This function, that can be wrapped within \code{\link{nosoiSim}}, runs a single-host transmission chain simulation, without any structure features in the host population. The simulation stops either at
#' the end of given time (specified by \code{length.sim}) or when the number of hosts infected threshold (\code{max.infected}) is crossed.
#'
#' @details
#' The \code{pExit} and \code{pTrans} functions should return a single probability (a number between 0 and 1), and \code{nContact} a positive natural number (positive integer) or 0.
#' @details
#' The \code{param} arguments should be a list of functions or NA. Each item name in the parameter list should have the same name as the argument in the corresponding function.
#' @details
#' The use of \code{timeDep} (switch to \code{TRUE}) makes the corresponding function use the argument \code{prestime} (for "present time").
#' @section Order of Arguments:
#' The user specified function's arguments should follow this order: \code{t} (mandatory), \code{prestime} (optional, only if timeDep is TRUE), \code{parameters} specified in the list.
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param nContact function that gives the number of potential transmission events per unit of time.
#' @param param.nContact parameter names (list of functions) for param.nContact.
#' @param timeDep.nContact is nContact dependent on the absolute time of the simulation? (TRUE/FALSE)
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pExit.
#' @param timeDep.pTrans is pTrans dependent on the absolute time of the simulation? (TRUE/FALSE)
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param timeDep.pExit is pExit dependent on the absolute time of the simulation? (TRUE/FALSE)
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param print.progress if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step print.progress is TRUE, step with which the progress message will be printed.
#'
#' @return An object of class \code{\link{nosoiSim}}, containing all results of the simulation.
#'
#' @seealso For simulations with a discrete structured host population, see \code{\link{singleDiscrete}}. For simulations with a structured population in continuous space, \code{\link{singleContinuous}}
#'
#' @examples
#' \donttest{
#'t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#'p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#'p_Exit_fct  <- function(t){return(0.08)}
#'
#'proba <- function(t,p_max,t_incub){
#'  if(t <= t_incub){p=0}
#'  if(t >= t_incub){p=p_max}
#'  return(p)
#'}
#'
#'time_contact <- function(t){round(rnorm(1, 3, 1), 0)}
#'
#'test.nosoi <- nosoiSim(type="single", popStructure="none",
#'                       length=40,
#'                       max.infected=100,
#'                       init.individuals=1,
#'                       nContact=time_contact,
#'                       param.nContact=NA,
#'                       pTrans = proba,
#'                       param.pTrans = list(p_max=p_max_fct,
#'                                           t_incub=t_incub_fct),
#'                       pExit=p_Exit_fct,
#'                       param.pExit=NA)
#'
#'test.nosoi
#'}
#' @export singleNone
singleNone <- function(length.sim,
                       max.infected,
                       init.individuals,
                       pExit, param.pExit, timeDep.pExit=FALSE,
                       nContact, param.nContact, timeDep.nContact=FALSE,
                       pTrans, param.pTrans, timeDep.pTrans=FALSE,
                       prefix.host="H",
                       print.progress=TRUE,
                       print.step=10){

  #Sanity check---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  #Parsing nContact
  nContactParsed <- parseFunction(nContact, param.nContact, as.character(quote(nContact)),timeDep=timeDep.nContact)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),timeDep=timeDep.pTrans)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),timeDep=timeDep.pExit)

  #Parsing all parameters
  ParamHost <- paramConstructor(param.pExit, param.pMove=NA, param.nContact, param.pTrans, param.sdMove=NA)

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(total.time = 1,
                             type = "single",
                             pop.A = nosoiSimOneConstructor(
                               N.infected = init.individuals,
                               table.hosts = iniTable(init.individuals, NA, prefix.host, ParamHost),
                               table.state = NA,
                               prefix.host = prefix.host,
                               popStructure = "none"))

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExitingMoving(res$host.info.A, pres.time, pExitParsed)

    res$host.info.A$table.hosts[exiting.full, `:=` (out.time = pres.time,
                                        active = FALSE)]

    if (!any(res$host.info.A$table.hosts[["active"]])) {break}

    #Step 1: Meeting & transmission ----------------------------------------------------

    df.meetTransmit <- meetTransmit(res$host.info.A, pres.time, positions = NULL, nContactParsed, pTransParsed)

    res$host.info.A <- writeInfected(df.meetTransmit, res$host.info.A, pres.time, ParamHost)

    if (print.progress == TRUE) progressMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time, print.step = print.step, length.sim = length.sim, max.infected.A = max.infected)
    if (res$host.info.A$N.infected > max.infected) {break}
  }

  endMessage(Host.count.A = res$host.info.A$N.infected, pres.time = pres.time)

  res$total.time <- pres.time

  return(res)
}
