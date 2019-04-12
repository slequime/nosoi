#' @title Parse function for later use
#'
#' @description
#' Parse a user-provided function to get the vectorized version.
#'
#' @param pFunc a function
#' @param param.pFunc a named list of arguments
#' @param name the name of the function
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

parseFunction <- function(pFunc, param.pFunc, name, diff=FALSE) {
  pFunc <- match.fun(pFunc)

  if ((diff == FALSE & length(formalArgs(pFunc)) > 1)|(diff == TRUE & length(formalArgs(pFunc)) > 2)) {
    if (!is.list(is.na(param.pFunc)) && is.na(param.pFunc)) {
      stop("There is a probleme with your function ", name, ": you should provide a parameter list named param.", name, ".")
    }

    if (is.list(param.pFunc)) {
      pFunc.param <- formalArgs(pFunc)[-1]
      if(! all(names(param.pFunc) %in% pFunc.param)) stop(paste0("Parameter name in param.", name, " should match the name used in ", name, "."))
    }
  }

    pFunc_eval <- function(prestime, inf.time,...) {
      t = prestime - inf.time
      x <- list(...)
      do.call(pFunc, c(list(t = t), x))
    }

    pFunc_vect <- function(prestime, parameters) {
      do.call(pFunc_eval, c(list(prestime = prestime), parameters))
    }

    pFunc_eval_args = c(formalArgs(pFunc_eval),formalArgs(pFunc)[-1])

    pFunc_eval_args = subset(pFunc_eval_args, pFunc_eval_args != "...")

    pFunc_vect_args = pFunc_eval_args[-1]

    return(list(vect = pFunc_vect,
                vectArgs = pFunc_vect_args))
  }
