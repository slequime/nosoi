#' @title get host info
#'
#' @description
#' Creates a \code{nosoiSim} object.
#'
#' @param res an object of class \code{nosoiSim}
#' @param what the information to get. One of "table.hosts", "N.infected", "table.state", "popStructure"
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with host informations
#'
#' @export
#'
##
getHostInfo <- function(res,
                        what = c("table.hosts", "N.infected", "table.state", "popStructure"),
                        pop = "A") {

  what <- match.arg(what)
  type <- res[["type"]]

  if(type == "single" && pop != "A") stop("There are no other hosts than 'A' in a single-host simulation.")

  if(what=="table.state" && pop == "A" && res$host.info.A$popStructure == "none") stop("There is no state information kept when the host population A has no structure.")
  if(what=="table.state" && pop == "B" && res$host.info.B$popStructure == "none") stop("There is no state information kept when the host population B has no structure.")

  if (pop == "A") return(res$host.info.A[[what]])
  if (pop == "B") return(res$host.info.B[[what]])

  stop(paste0("Population ", pop, " is not recognized."))
}

#' @title Get table hosts
#'
#' @description
#' Get "table.hosts"
#' TODO: describe the table
#'
#' @param res an object of class \code{nosoiSim}
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with hosts table informations
#'
##
getTableHosts <- function(res, pop = "A") {
  return(getHostInfo(res, "table.hosts", pop))
}

#' @title Get table states
#'
#' @description
#' Get "table.state"
#' TODO: describe the table
#'
#' @param res an object of class \code{nosoiSim}
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return a data.table with state table informations
#'
##
getTableState <- function(res, pop = "A") {
  return(getHostInfo(res, "table.state", pop))
}

