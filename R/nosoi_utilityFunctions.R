
#' @title Parse function for later use
#'
#' @description
#' Parse a user-provided function to get the vectorized version.
#'
#' @param pFunc a function
#' @param param.pFunc a named list of arguments
#' @param name the name of the function
#' @param diff is the function differential according to state/env.variable? (TRUE/FALSE)
#' @param timeDep is the function differential according to absolute time? (TRUE/FALSE)
#' @param hostCount is the function differential according to host count? (TRUE/FALSE)
#' @param continuous is the function to be used in a continuous space? (TRUE/FALSE)
#' @param stateNames name of the states (vector) in case of discrete structure.
#'
#' @return list of parsed quantities:
#' \itemize{
#'  \item{"vect"}{Vectorized version of the function.}
#'  \item{"vectArgs"}{Vector of arguments for the vectorized function.}
#' }
#'
#' @keywords internal
##

parseFunction <- function(pFunc, param.pFunc, name, diff=FALSE, timeDep=FALSE, hostCount=FALSE, continuous=FALSE, stateNames=NA) {

  FunctionSanityChecks(pFunc, name, param.pFunc, timeDep, diff, hostCount, continuous, stateNames)

  pFunc <- match.fun(pFunc)

  if (timeDep == FALSE){
    pFunc_eval <- function(prestime, inf.time,...) {
      t = prestime - inf.time
      x <- list(...)
      do.call(pFunc, c(list(t = t), x))
    }

    pFunc_eval_args = c(formalArgs(pFunc_eval),formalArgs(pFunc)[-1])
  }

  if (timeDep == TRUE){
    pFunc_eval <- function(prestime, inf.time,...) {
      t = prestime - inf.time
      x <- list(...)
      do.call(pFunc, c(list(t = t),list(prestime=prestime), x))

    }
    pFunc_eval_args = c(formalArgs(pFunc_eval),formalArgs(pFunc)[c(-1,-2)])
  }

  pFunc_eval_args = subset(pFunc_eval_args, pFunc_eval_args != "...")

  pFunc_vect_args = pFunc_eval_args[-1]

  pFunc_vect <- function(prestime, parameters) {
    do.call(pFunc_eval, c(list(prestime = prestime), parameters))
  }

  return(list(vect = pFunc_vect,
              vectArgs = pFunc_vect_args))
}

## Draw K bernouillis with probas p1, ..., pK

#' @title Draw newly infected
#'
#' @description
#' For each encounter, simulate wether a new individual is infected.
#'
#' @param p vector of size K, giving the probability that each encounter
#' leads to an infection.
#'
#' @return Boolean vector giving the newly infected individuals.
#' @keywords internal
##

drawBernouilli <- function(p) {
  if (!is.vector(p)) stop("Function 'drawBernouilli' should be applied to a vector.")
  return(runif(length(p), 0, 1) < p)
}

#' @title get Position Infected
#'
#' @description
#' Return the relevent position of the infected individual.
#'
#' @param nosoiSim an object of class \code{nosoiSim}.
#' @param df.meetTransmit current data.table for transmission
#' @param i individual
#'
#' @return A vector of positions
#'
#' @keywords internal
##
getPositionInfected <- function(nosoiSim, df.meetTransmit, i) {
  if (nosoiSim$popStructure == "none") return(NA)
  if (nosoiSim$popStructure == "discrete") return(df.meetTransmit[i, ]$current.in)
  if (nosoiSim$popStructure == "continuous") return(c(df.meetTransmit[i, ]$current.in.x, df.meetTransmit[i, ]$current.in.y))
  stop(paste0("Geographical structure ", nosoiSim$popStructure, " is not implemented."))
}

#' @title Should we build the table.host table
#'
#' @param nosoiSim an object of class \code{nosoiSim}.
#'
#' @return boolean
#'
#' @keywords internal
##
keepState <- function(nosoiSim) {
  if (nosoiSim$popStructure == "none") return(FALSE)
  if (nosoiSim$popStructure == "discrete") return(TRUE)
  if (nosoiSim$popStructure == "continuous") return(TRUE)
  stop(paste0("Geographical structure \"", nosoiSim$popStructure, "\" is not implemented."))
}

#' @title Param concatenator
#'
#' @description
#' Creates a \code{ParamHost} object.
#'
#' @inheritParams singleContinuous
#'
#' @return An object of class \code{ParamHost}
#'
#' @keywords internal
##

paramConstructor <- function(param.pExit, param.pMove, param.nContact, param.pTrans,
                             param.sdMove) {

  checkna <- function(param){return((length(param) == 1) && is.na(param))}

  if (checkna(param.pExit)) param.pExit <- NULL
  if (checkna(param.pMove)) param.pMove <- NULL
  if (checkna(param.nContact)) param.nContact <- NULL
  if (checkna(param.sdMove)) param.sdMove <- NULL
  if (checkna(param.pTrans)) param.pTrans <- NULL

  merged <- c(param.pTrans, param.pMove, param.nContact, param.pExit, param.sdMove)

  ParamHost <- merged[!duplicated(merged)]

  if(!is.null(ParamHost)) class(ParamHost) <- "ParamHost"

  return(ParamHost)
}

#' @title Apply a function to table.host
#'
#' @description
#' Return a vector of the results of the function
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @param pres.time current time
#' @param pasedFunction parsed exit/moving function
#' @param active.hosts a boolean vector of active hosts
#'
#' @return result vector
#'
#' @keywords internal
##
applyFunctionToHosts <- function(res, pres.time, pasedFunction, active.hosts) {
  return(res$table.hosts[active.hosts,
                         pasedFunction$vect(prestime = pres.time, .SD),
                         by="hosts.ID",
                         .SDcols = pasedFunction$vectArgs][["V1"]])
}

#' @title Get Exiting or Moving individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @param pres.time current time
#' @param pasedFunction parsed exit/moving function
#'
#' @return Boolean vector
#'
#' @keywords internal
##
getExitingMoving <- function(res, pres.time, pasedFunction) {

  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

  if (any(active.hosts)) {

    p.exitMove.values <- applyFunctionToHosts(res, pres.time, pasedFunction, active.hosts)

    exitMove <- drawBernouilli(p.exitMove.values) #Draws K bernouillis with various probability (see function for more detail)
  }

  if(all(active.hosts == FALSE)) exitMove <- FALSE

  exitMove.full <- active.hosts
  exitMove.full[exitMove.full] <- exitMove

  return(exitMove.full)
}

#' @title Summarise position of hosts in a discrete or discretized (raster) space
#'
#' @description
#' Returns a table of the number of hosts per discrete or discretized (raster cells) states.
#'
#' @param res an object of class \code{nosoiSimOne}.
#'
#' @return data.table containing state.ID (col1) and count (col2)
#'
#' @keywords internal
##

updateHostCount <- function(res, res.B=NULL, type) {

  #To avoid notes (use of data.table shortcuts)
  count.A <- NULL
  count.B <- NULL
  current.cell.raster <- NULL
  current.in <- NULL

  #New IDEA: create count table for both A & B at same time with both needed states => merge

  #Step 1: count Table
  active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)
  if(type == "discrete"){host.table.count.A <- as.data.table(table(res$table.hosts[active.hosts]$current.in))}
  if(type == "continuous"){host.table.count.A <- as.data.table(table(res$table.hosts[active.hosts]$current.cell.raster))}
  host.table.count.A$V1 = as.character(host.table.count.A$V1)
  colnames(host.table.count.A) = c("state.ID","count")

  if(!is.null(res.B)){
    active.hosts.B <- res.B$table.hosts[["active"]] == 1 #active hosts (boolean vector)

    if(type == "discrete"){host.table.count.B <- as.data.table(table(res.B$table.hosts[active.hosts.B]$current.in))}
    if(type == "continuous"){host.table.count.B <- as.data.table(table(res.B$table.hosts[active.hosts.B]$current.cell.raster))}
    host.table.count.B$V1 = as.character(host.table.count.B$V1)
    colnames(host.table.count.B) = c("state.ID","count")
  }

  if(is.null(res.B)) {host.table.count.B = data.table(V1=0,N=0)[-1]
  host.table.count.B$V1 = as.character(host.table.count.B$V1)
  colnames(host.table.count.B) = c("state.ID","count")}

  host.table.count.tot = merge(host.table.count.A,host.table.count.B,by="state.ID", all=TRUE, suffixes = c(".A",".B"))
  setkey(host.table.count.tot,"state.ID")

  #Change NA into 0
  host.table.count.tot[is.na(count.B), `:=` (count.B = 0)]
  host.table.count.tot[is.na(count.A), `:=` (count.A = 0)]

  #Step 2: update counts

  if(type == "discrete"){
    if(is.null(res.B)){
      for (i in host.table.count.tot$state.ID){
        res$table.hosts[(active.hosts & current.in == i), `:=` (host.count = host.table.count.tot[i][["count.A"]])]
      }
    }

    if(!is.null(res.B)){
      for (i in host.table.count.tot$state.ID){

        res$table.hosts[(active.hosts & current.in == i), `:=` (host.count = host.table.count.tot[i][["count.A"]],
                                                                host.count.B = host.table.count.tot[i][["count.B"]])]

        res.B$table.hosts[(active.hosts.B & current.in == i), `:=` (host.count = host.table.count.tot[i][["count.A"]],
                                                                    host.count.B = host.table.count.tot[i][["count.B"]])]
      }
    }
  }

  if(type == "continuous"){
    if(is.null(res.B)){
      for (i in host.table.count.tot$state.ID){
        res$table.hosts[(active.hosts & current.cell.raster == i), `:=` (host.count = host.table.count.tot[as.character(i)][["count.A"]])]
      }
    }

    if(!is.null(res.B)){
      for (i in host.table.count.tot$state.ID){
        res$table.hosts[(active.hosts & current.cell.raster == i), `:=` (host.count = host.table.count.tot[as.character(i)][["count.A"]],
                                                                host.count.B = host.table.count.tot[as.character(i)][["count.B"]])]

        res.B$table.hosts[(active.hosts.B & current.cell.raster == i), `:=` (host.count = host.table.count.tot[as.character(i)][["count.A"]],
                                                                             host.count.B = host.table.count.tot[as.character(i)][["count.B"]])]
      }
    }
  }
}
