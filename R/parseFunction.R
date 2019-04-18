#' @title Parse function for later use
#'
#' @description
#' Parse a user-provided function to get the vectorized version.
#'
#' @param pFunc a function
#' @param param.pFunc a named list of arguments
#' @param name the name of the function
#' @param diff is the function differential according to state/env.variable? (TRUE/FALSE)
#' @param timeDep
#' @param continuous
#' @param stateNames
#'
#' @return list of parsed quantities:
#' \itemize{
#'  \item{"type"}{A string, "simple" of "complex": are there several arguments to the function ?}
#'  \item{"nArgs"}{Number of arguments (besides time).}
#'  \item{"vect"}{Vectorized version of the function.}
#'  \item{"vectArgs"}{Vector of arguments for the vectorized function.}
#' }
#'
#' @keywords internal
##

parseFunction <- function(pFunc, param.pFunc, name, diff=FALSE, timeDep=FALSE, continuous=FALSE, stateNames=NA) {

  FunctionSanityChecks(pFunc, name, param.pFunc, timeDep, diff, continuous, stateNames)

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
