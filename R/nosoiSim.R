nosoiSim <- function(length=NA, max.infected=NA,type="single",geo="none",parallel=FALSE,...){
library(data.table)
  #Sanity checks -------------
  if (is.na(length) | length <= 1) message("You must specify a length (in time units) for your simulation.")
  if (is.na(max.infected) | max.infected <= 1) message("You must specify a maximum number of infected hosts.")
  if (! type %in% c("single","dual")) message("Type of transmission should be 'single' or 'dual'-host.")
  if (! geo %in% c("none","discrete","continuous")) message("Unrecognized parameters for location-explicit model")
  if (! parallel %in% c(TRUE,FALSE)) message("Parallel parameter should be TRUE or FALSE")

  #Loading correct script ------------------
  if(type=="single" & geo=="none" & parallel==FALSE) {
    source("./single-none-nonparallel.R")
    single_none_parallel(length=length,max.infected=max.infected,...)
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
}
