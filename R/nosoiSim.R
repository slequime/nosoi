#' Top-level function to use Nosoi
#'
#' @param type specifies which type of pathogen we are interested in, either "single" or "dual"-host (e.g. arboviruses).
#' @param structure specifies if the population in which the transmission is to occur is structured (discrete variable).
#' @param continuous specifies if the spread is to occur in a continuous space (continuous diffusion).
#' @param ... arguments to be passed on to the simulator (see below).
#'
#' @details This function determines which general settings the user wants to use for his simulation.
#' @details All other arguments are passed down to the chosen simulator itself, such as \code{\link{singleNone}}
#' @export nosoiSim
#' @import data.table
#' @import utils
#' @import methods
#' @import stats

nosoiSim <- function(type="single", structure=FALSE, continuous=FALSE, ...){
  #Sanity checks -------------
  if (! type %in% c("single","dual")) stop("Type of transmission should be 'single' or 'dual'-host.")
  if (! structure %in% c(TRUE,FALSE)) stop("Unrecognized parameters for population structure, should be TRUE or FALSE.")

  #Loading correct script ------------------
  if (type=="single" && structure==FALSE) {
    output <- singleNone(...)
  }
  if (type=="single" && structure==TRUE && continuous == FALSE) {
    output <- singleDiscrete(...)
  }
  if (type=="single" && structure==TRUE && continuous == TRUE) {
    output <- singleContinuous(...)
  }

  if (type=="dual" && structure==FALSE) {
    output <- dualNone(...)
  }
  if (type=="dual" && structure==TRUE && continuous == FALSE) {
    output <- dualDiscrete(...)
  }

  #To be implemented
  if(type=="dual" && structure==TRUE && continuous == TRUE) stop("This version has not been implemented yet.")

  return(output)
}
