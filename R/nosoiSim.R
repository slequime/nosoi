#' Top-level function to use Nosoi
#'
#' @param type specifies which type of pathogen we are interested in, either "single" or "dual"-host (e.g. arboviruses).
#' @param structure specifies if the population in which the transmission is to occur is structured (discrete variable).
#' @param ... arguments to be passed on to the simulator (see below).
#'
#' @details This function determines which general settings the user wants to use for his simulation.
#' @details All other arguments are passed down to the chosen simulator itself, such as \code{\link{SingleNone}}
#' @export nosoiSim
#' @import data.table
#' @import utils
#' @import methods
#' @import stats

nosoiSim <- function(type="single",structure=FALSE,...){
  #Sanity checks -------------
  if (! type %in% c("single","dual")) message("Type of transmission should be 'single' or 'dual'-host.")
  if (! structure %in% c(TRUE,FALSE)) message("Unrecognized parameters for population structure, should be TRUE or FALSE.")

  #Loading correct script ------------------
  if(type=="single" & structure==FALSE) {
    output = SingleNone(...)
  }
  #To be implemented
  if(type=="single" & structure==TRUE) message("This version has not been implemented yet.")
  if(type=="dual" & structure==FALSE) message("This version has not been implemented yet.")
  if(type=="dual" & structure==TRUE) message("This version has not been implemented yet.")
  return(output)
}
