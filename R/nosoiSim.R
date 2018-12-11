#' Top-level function to use Nosoi
#'
#' @param type specifies which type of pathogen we are interested in, either "single" or "dual"-host (e.g. arboviruses).
#' @param geo specifies if Nosoi should be spatially explicit or not ("none"), and if yes, which type of geographic space ("continuous" or "discrete").
#' @param parallel if simulation should be run in parallel or not.
#' @param ... arguments to be passed on to the simulator (see below).
#'
#' @details This function determines which general settings the user wants to use for his simulation.
#' @details All other arguments are passed down to the chosen simulator itself, such as \code{\link{single_none_parallel}}
#' @export nosoiSim
#' @import data.table

nosoiSim <- function(type="single",geo="none",parallel=FALSE,...){
  #Sanity checks -------------
  if (! type %in% c("single","dual")) message("Type of transmission should be 'single' or 'dual'-host.")
  if (! geo %in% c("none","discrete","continuous")) message("Unrecognized parameters for location-explicit model")
  if (! parallel %in% c(TRUE,FALSE)) message("Parallel parameter should be TRUE or FALSE")

  #Loading correct script ------------------
  if(type=="single" & geo=="none" & parallel==FALSE) {
    output = single_none_parallel(length=length,max.infected=max.infected,...)
}
  #To be implemented
  if(type=="single" & geo=="none" & parallel==TRUE) message("This version has not been implemented yet.")
  if(type=="single" & geo=="discrete" & parallel==FALSE) message("This version has not been implemented yet.")
  if(type=="single" & geo=="discrete" & parallel==TRUE) message("This version has not been implemented yet.")
  if(type=="single" & geo=="continuous" & parallel==FALSE) message("This version has not been implemented yet.")
  if(type=="single" & geo=="continuous" & parallel==TRUE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="none" & parallel==FALSE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="none" & parallel==TRUE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="discrete" & parallel==FALSE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="discrete" & parallel==TRUE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="continuous" & parallel==FALSE) message("This version has not been implemented yet.")
  if(type=="dual" & geo=="continuous" & parallel==TRUE) message("This version has not been implemented yet.")
return(output)
  }
