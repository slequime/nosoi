#This file holds all the function directly related to the nosoiSim function, including the object constructors

#' @title Top-level function to use \code{nosoi}.
#'
#' @description This function determines which general settings the user wants to use for his simulation. All other arguments are passed down to the chosen simulator itself, such as \code{\link{singleNone}}, \code{\link{singleDiscrete}}, \code{\link{singleContinuous}}, \code{\link{dualNone}}, \code{\link{dualDiscrete}} or \code{\link{dualContinuous}}.
#'
#' @param type specifies which type of pathogen we are interested in, either "single" or "dual"-host (e.g. arboviruses).
#' @param popStructure specifies if the population in which the transmission is to occur is structured ("none", "discrete" or "continuous").
#' @param ... arguments to be passed on to the chosen simulator itself, such as \code{\link{singleNone}}, \code{\link{singleDiscrete}}, \code{\link{singleContinuous}}, \code{\link{dualNone}}, \code{\link{dualDiscrete}} or \code{\link{dualContinuous}}.
#'
#' @return An object of class \code{nosoiSim}, containing all results of the simulation. Class \code{nosoiSim} object have the following slots:
#' \describe{
#'    \item{total.time}{Number of time steps the simulation ran (integer).}
#'
#'    \item{type}{String giving the simulation type ("single" or "dual" host).}
#'
#'    \item{host.info.A: object of class \code{nosoiSimOne}}{
#'        \describe{
#'           \item{N.infected}{Number of infected hosts (integer).}
#'           \item{table.hosts}{Table containing the results of the simulation  (see \code{\link{getTableHosts}} for more details on the table).}
#'           \item{table.state}{Table containing the results of the simulation, focusing on the movement history of each host (see \code{\link{getTableState}} for more details on the table).}
#'           \item{prefix.host}{String containing the prefix used to name hosts (character string).}
#'           \item{popStructure}{String giving the population structure (one of "none", "discrete" or "continuous").}
#'    }}
#'
#'    \item{host.info.B: object of class \code{nosoiSimOne}}{
#'    Same structure as \code{host.info.A}, but for host B (if it exists).
#'    }
#'    }
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
#'time_contact = function(t){round(rnorm(1, 3, 1), 0)}
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
#'test.nosoi
#'}
#'
#' @seealso
#' \describe{
#'   \item{Individual simulation functions:}{
#'     \code{\link{singleNone}}, \code{\link{singleDiscrete}},
#'     \code{\link{singleContinuous}}, \code{\link{dualNone}},
#'     \code{\link{dualDiscrete}} and \code{\link{dualContinuous}}.
#'     }
#'     \item{Functions to extract the results:}{
#'     \code{\link{getTableHosts}}, \code{\link{getTableState}}
#'     }
#'     \item{Summary statistics functions:}{
#'     \code{\link{nosoiSummary}},
#'     \code{\link{getCumulative}}, \code{\link{getDynamic}}, \code{\link{getR0}}
#'     }
#' }
#'
#' @export nosoiSim
#'
#' @import data.table
#' @import methods
#' @import stats

nosoiSim <- function(type="single", popStructure="none", ...){

  #Sanity checks -------------
  if (! type %in% c("single","dual")) stop("Type of transmission should be 'single' or 'dual'-host.")
  if (! popStructure %in% c("none","discrete","continuous")) stop("Unrecognized parameters for population structure, should be 'none,'discrete' or 'continuous'.")

  #Loading correct script ------------------
  if (type=="single" && popStructure=="none") {
    output <- singleNone(...)
  }
  if (type=="single" && popStructure=="discrete") {
    output <- singleDiscrete(...)
  }
  if (type=="single" && popStructure=="continuous") {
    output <- singleContinuous(...)
  }

  if (type=="dual" && popStructure=="none") {
    output <- dualNone(...)
  }
  if (type=="dual" && popStructure=="discrete") {
    output <- dualDiscrete(...)
  }

  if(type=="dual" && popStructure=="continuous") {
    output <- dualContinuous(...)
  }

  return(output)
}

#' @title nosoiSim Constructor
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param pres.time current time of the simulation
#' @param type population structure (one of "single or "dual)
#' @param pop.A an object of class \code{nosoiSimOne} for population A
#' @param pop.B an object of class \code{nosoiSimOne} for population B
#'
#' @return An object of class \code{nosoiSim}
#'
#' @keywords internal
##

nosoiSimConstructor <- function(total.time,
                                type = c("single", "dual"),
                                pop.A,
                                pop.B = NULL) {

  type <- match.arg(type)

  res <- list(total.time = total.time,
              type = type,
              host.info.A = pop.A,
              host.info.B = pop.B)

  class(res) <- "nosoiSim"

  return(res)

}

##
#' @export
#' @method print nosoiSim
##
print.nosoiSim <- function(x, ...){
  cat("A nosoiSim object, representing a simulated epidemy ")
  cat(paste0("for a ", x$type, " host with "))
  if (x$host.info.A$popStructure == "none") cat("no structure.")
  if (x$host.info.A$popStructure == "discrete") cat("a discrete structure.")
  if (x$host.info.A$popStructure == "continuous") cat("a continuous structure.")
  cat(endMessageText(x$host.info.A$N.infected, x$host.info.A$N.infected,
                     x$total.time, x$type, ""))
  cat("\nUse function 'summary' for summary statistics, and functions 'getTableHosts' and 'getTableState' to extract the generated data.")
}

#' @title nosoiSimOne Constructor
#'
#' @description
#' Creates a \code{nosoiSimOne} object.
#'
#' @param N.infected number of infected hosts
#' @param table.hosts data.table of hosts
#' @param table.state data.table of hosts movement
#' @param popStructure geographical structure (one of "none, "discrete" or "continuous")
#'
#' @return An object of class \code{nosoiSimOne}
#'
#' @keywords internal
##
nosoiSimOneConstructor <- function(N.infected, table.hosts, table.state, prefix.host,
                                   popStructure = c("none", "discrete", "continuous")) {

  popStructure <- match.arg(popStructure)

  res <- list(N.infected = N.infected,
              table.hosts = table.hosts,
              table.state = table.state,
              prefix.host = prefix.host,
              popStructure = popStructure)

  class(res) <- "nosoiSimOne"

  return(res)

}
