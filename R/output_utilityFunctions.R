#' @title Extracts specific data from a \code{nosoiSim} object
#'
#' @description
#' This function extracts data user-defined data (i.e. \code{table.hosts}, \code{N.infected}, \code{table.state} or \code{popStructure}) from a \code{\link{nosoiSim}} object.
#'
#' @param nosoi.output an object of class \code{nosoiSim}
#' @param what the data to get, among \code{table.hosts}, \code{N.infected}, \code{table.state} or \code{popStructure}.
#' @param pop the population to be extracted (one of "A" or "B")
#'
#' @return Returns a \code{\link[data.table:data.table-package]{data.table}} with the requested data.
#'
#' @seealso To directly extract \code{table.hosts} or \code{table.state}, you can also use \code{\link{getTableHosts}} and \code{\link{getTableState}} respectively.
#' @examples
#' \donttest{
#' t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#' p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#' p_Exit_fct  <- function(t){return(0.08)}
#'
#' proba <- function(t,p_max,t_incub){
#'   if(t <= t_incub){p=0}
#'   if(t >= t_incub){p=p_max}
#'   return(p)
#' }
#'
#' time_contact <- function(t){round(rnorm(1, 3, 1), 0)}
#'
#' test.nosoi <- nosoiSim(type="single", popStructure="none",
#'                        length=40,
#'                        max.infected=100,
#'                        init.individuals=1,
#'                        nContact=time_contact,
#'                        param.nContact=NA,
#'                        pTrans = proba,
#'                        param.pTrans = list(p_max=p_max_fct,
#'                                            t_incub=t_incub_fct),
#'                        pExit=p_Exit_fct,
#'                        param.pExit=NA)
#'
#'
#' data.extracted <- getHostData(test.nosoi, "table.hosts", "A")
#'}
#' @export getHostData

getHostData <- function(nosoi.output,
                        what = c("table.hosts", "N.infected", "table.state", "popStructure"),
                        pop = "A") {

  what <- match.arg(what)
  type <- nosoi.output[["type"]]

  if(type == "single" && pop != "A") stop("There are no other hosts than 'A' in a single-host simulation.")

  if(what=="table.state" && pop == "A" && nosoi.output$host.info.A$popStructure == "none") stop("There is no state information kept when the host population A has no structure.")
  if(what=="table.state" && pop == "B" && nosoi.output$host.info.B$popStructure == "none") stop("There is no state information kept when the host population B has no structure.")

  if (pop == "A") return(nosoi.output$host.info.A[[what]])
  if (pop == "B") return(nosoi.output$host.info.B[[what]])

  stop(paste0("Population ", pop, " is not recognized."))
}

#' @title Extracts \code{table.hosts} from a \code{nosoiSim} object
#'
#' @description This function extracts the \code{table.hosts} for the request host-type from a \code{\link{nosoiSim}} object.
#'
#' @param nosoi.output an object of class \code{\link{nosoiSim}}
#' @param pop the host-type to be extracted (either "A" or "B", if not dual-host, then "A")
#'
#' @return Returns a \code{\link[data.table:data.table-package]{data.table}} with the requested data. The \code{table.hosts} (class \code{\link[data.table:data.table-package]{data.table}}) contains informations about each host that has been simulated (one row is one host).
#' The structure of the table is the following:
#' \describe{
#'    \item{hosts.ID}{Unique identifier for the host, based on user-defined prefix and an integer.}
#'    \item{inf.by}{Unique identifier for the host that infected the current one.}
#'    \item{inf.in}{(only if structure is present) State or coordinates (in that case inf.in.x and inf.in.y) in which the host was infected.}
#'    \item{current.in}{(only if structure is present) State or coordinates (in that case current.in.x and current.in.y) in which the host is at the end of the simulation.}
#'    \item{current.env.value}{(only if continuous structure is present) Environmental value (raster cell value) in which the host is at the end of the simulation.}
#'    \item{current.cell.raster}{(only if continuous structure is present) Raster cell numeric ID in which the host is at the end of the simulation.}
#'    \item{host.count}{(only if structure is present) Host count in the current state or raster cell (beware, updated only if used).}
#'    \item{inf.time}{When did the host enter the simulation (infection time).}
#'    \item{out.time}{When did the host exit the simulation (NA if still active).}
#'    \item{active}{Is the host still active at the end of the simulation (TRUE for YES, FALSE for NO)?}
#'    \item{parameters}{The remaining columns are the sampled values for the individual-based parameters (if any) specified by the user.}
#'    }
#'
#'
#' @export getTableHosts
getTableHosts <- function(nosoi.output, pop = "A") {
  return(getHostData(nosoi.output, "table.hosts", pop))
}

#' @title Extracts \code{table.state} from a \code{nosoiSim} object
#'
#' @description This function extracts the \code{table.state} for the request host-type from a \code{\link{nosoiSim}} object. \code{table.state} is present only if there is any structure (discrete or continuous) used.
#'
#' @param nosoi.output an object of class \code{\link{nosoiSim}}
#' @param pop the host-type to be extracted (either "A" or "B", if not dual-host, then "A")
#'
#' @return Returns a \code{\link[data.table:data.table-package]{data.table}} with the requested data. The \code{table.state} (class \code{\link[data.table:data.table-package]{data.table}}) contains informations the location of each host during time (one row is one host at one location).
#' The structure of the table is the following:
#' \describe{
#'    \item{hosts.ID}{Unique identifier for the host, based on user-defined prefix and an integer.}
#'    \item{state}{State or coordinates (in that case state.x and state.y) in which the host is during that time interval.}
#'    \item{current.env.value}{(only if continuous structure is present) Environmental value (raster cell value) in which the host is at the end of the simulation.}
#'    \item{current.cell.raster}{(only if continuous structure is present) Raster cell numeric ID in which the host is at the end of the simulation.}
#'    \item{time.from}{Time-step at which the host moved to the location.}
#'    \item{time.to}{Time-step at which the host exited the location (either by exiting the simulation or moving somewhere else).}
#'    }
#'
#' @export getTableState

getTableState <- function(nosoi.output, pop = "A") {
  return(getHostData(nosoi.output, "table.state", pop))
}

