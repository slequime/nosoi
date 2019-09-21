
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

  if (hostCount && any("host.count" == formalArgs(pFunc))) { ## Replace user defined "host.count" by "host.count.A"
    ff <- formals(pFunc)
    names(ff)["host.count" == names(ff)] <- "host.count.A"
    formals(pFunc) <- ff
    body(pFunc) <- substituteDirect(body(pFunc), list(host.count = substitute(host.count.A)))
  }

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
#' For each encounter, simulate whether a new individual is infected.
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
#' Return the relevant position of the infected individual.
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
  if (nosoiSim$popStructure == "discrete") return(df.meetTransmit[[i, "current.in"]])
  if (nosoiSim$popStructure == "continuous") return(c(df.meetTransmit[[i, "current.in.x"]], df.meetTransmit[[i, "current.in.y"]]))
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
                         by = "hosts.ID",
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

  active.hosts <- res$table.hosts[["active"]] #active hosts (boolean vector)

  if (any(active.hosts)) {

    p.exitMove.values <- applyFunctionToHosts(res, pres.time, pasedFunction, active.hosts)

    exitMove <- drawBernouilli(p.exitMove.values) #Draws K bernouillis with various probability (see function for more detail)
  }

  if (!any(active.hosts)) exitMove <- FALSE

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

  if(type == "discrete"){
    name_cur <- "current.in"
  } else if (type == "continuous") {
    name_cur <- "current.cell.raster"
  }

  active.hosts <- res$table.hosts[["active"]] #active hosts (boolean vector)

  res$table.hosts[active.hosts, by=name_cur, `:=` (host.count.A = .N)]

  if(!is.null(res.B)){
    active.hosts.B <- res.B$table.hosts[["active"]] #active hosts (boolean vector)

    res.B$table.hosts[active.hosts.B, by=name_cur, `:=` (host.count.B = .N)]

    res$table.hosts[active.hosts, by=name_cur,
                    `:=` (host.count.B = get_other_count(table_hosts = res.B$table.hosts,
                                                         name_cur = name_cur,
                                                         active_hosts = active.hosts.B,
                                                         name_count = "host.count.B",
                                                         .SD)),
                    .SDcols = name_cur]
    res.B$table.hosts[active.hosts.B, by=name_cur,
                      `:=` (host.count.A = get_other_count(table_hosts = res$table.hosts,
                                                         name_cur = name_cur,
                                                         active_hosts = active.hosts,
                                                         name_count = "host.count.A",
                                                         .SD)),
                      .SDcols = name_cur]
  }
}

#' @title Get host count from table
#'
#' @description
#' Returns the host count for a table and a state
#'
#' @param table_hosts a table.hosts
#' @param name_cur name of the cell, one of "current.in" or "current.cell.raster"
#' @param active_hosts a boolean vector of active hosts
#' @param name_count name of the count, one of "host.count.A" or "host.count.B"
#' @param cur a data.table with the state
#'
#' @return data.table containing state.ID (col1) and count (col2)
#'
#' @keywords internal
##
get_other_count <- function(table_hosts, name_cur, active_hosts, name_count, cur) {
  if (nrow(cur) == 0) return(0L)
  val <- table_hosts[active_hosts & get(name_cur) == cur[[1,1]], name_count, with=FALSE]
  return(ifelse(nrow(val) == 0, 0L, val[[1, 1]]))
}
