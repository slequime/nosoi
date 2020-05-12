#This file holds all the function directly related to SanityChecks and messages to the user.

#' @title Checks if the simulator can start
#'
#' @description
#' Checks if the simulator can start: did the user provide a length, a maximum number of infected individuals and a number of "seeding" individuals?
#'
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals  number of initially infected individuals.
#'
#' @keywords internal
##

CoreSanityChecks <- function(length.sim, max.infected, init.individuals) {
  if (is.na(length.sim) || length.sim <= 1) stop("You must specify a length (in time units) for your simulation (bigger than 1).")
  if (is.na(max.infected) || max.infected <= 1) stop("You must specify a maximum number of infected hosts (bigger than 1).")
  if (is.na(init.individuals) || init.individuals < 1 || !init.individuals%%1==0) stop("The transmission chain should be started by 1 or more (integer) individuals.")
}

#' @title Checks if a function is properly formatted
#'
#' @description
#' Checks if the function is properly formatted given the user's input and nosoi's requirements.
#'
#' @param pFunc a function
#' @param name the name of the function
#' @param param.pFunc a named list of arguments
#' @param timeDep is the function differential according to absolute time? (TRUE/FALSE)
#' @param diff is the function differential according to state/env.variable? (TRUE/FALSE)
#' @param hostCount is the function differential according to host count? (TRUE/FALSE)
#' @param structure is the function to be used in a structured population? (TRUE/FALSE)
#' @param continuous is the function to be used in a continuous space? (TRUE/FALSE)
#' @param stateNames name of the states (vector) in case of discrete structure.
#'
#' @keywords internal
##

FunctionSanityChecks <- function(pFunc, name, param.pFunc, timeDep, diff, hostCount, continuous, stateNames) {

  if (!is.function(pFunc)) stop("You must specify ",name, " as a function.")

  pFunc <- match.fun(pFunc)
  if (!formalArgs(pFunc)[1] == "t") stop(name, " must be a function of 't', placed as a first argument of the function.")

  if(hostCount == TRUE && diff == FALSE) stop("diff.", name, " should be TRUE to use hostCount.", name, ".")

  if ((diff == FALSE && timeDep == FALSE && length(formalArgs(pFunc)) > 1)
      || (diff == TRUE && timeDep == FALSE && hostCount == FALSE && length(formalArgs(pFunc)) > 2)
      || (diff == TRUE && timeDep == FALSE && hostCount == TRUE && length(formalArgs(pFunc)) > 3)
      || (diff == FALSE && timeDep == TRUE && length(formalArgs(pFunc)) > 2)
      || (diff == TRUE && timeDep==TRUE && hostCount == FALSE && length(formalArgs(pFunc)) > 3)
      || (diff == TRUE && timeDep==TRUE && hostCount == TRUE && length(formalArgs(pFunc)) > 4)) {
    if (!is.list(param.pFunc) && is.na(param.pFunc)) {
      stop("There is a probleme with your function ", name, ": you should provide a parameter list named param.", name, ".")
    }

    if (is.list(param.pFunc)) {
      pFunc.param <- formalArgs(pFunc)[-1]
      if(! all(names(param.pFunc) %in% pFunc.param)) stop("Parameter name in param.", name, " should match the name used in ", name, ".")
    }
  }

  if((diff == TRUE || timeDep == TRUE) && length(formalArgs(pFunc)) < 2) stop("Your are missing some function argument in ",name,". diff and/or timeDep.", name, " is/are TRUE.")

  if (diff == TRUE && continuous == FALSE){
    if (any(!sapply(paste0('current.in == "',stateNames,'"'),
                   function(pat) grepl(pat, paste0(as.character(body(pFunc)), collapse = " "))))) stop(name, " should have a realisation for each possible state. diff.", name, " is TRUE.")
  }

  if(timeDep == TRUE && formalArgs(pFunc)[2] != "prestime") stop(name, " should have 'prestime' as the second variable. timeDep.", name, " is TRUE.")

  if(timeDep == FALSE && diff == TRUE && continuous == FALSE && formalArgs(pFunc)[2] != "current.in") stop(name, " should have 'current.in' as the second variable. diff.", name, " is TRUE.")
  if(timeDep == TRUE && diff == TRUE && continuous == FALSE && formalArgs(pFunc)[3] != "current.in") stop(name, " should have 'current.in' as the third variable. diff.", name, " is TRUE.")

  if(timeDep == FALSE && diff == TRUE && continuous == TRUE && formalArgs(pFunc)[2] != "current.env.value") stop(name, " should have 'current.env.value' as the second variable. diff.", name, " is TRUE.")
  if(timeDep == TRUE && diff == TRUE && continuous == TRUE && formalArgs(pFunc)[3] != "current.env.value") stop(name, " should have 'current.env.value' as the third variable. diff.", name, " is TRUE.")
}

#' @title Checks if the matrix is properly formatted
#'
#' @description
#' Checks if the transition matrix is properly formatted given the user's input and nosoi's requirements.
#'
#' @param structure.matrix transition matrix (probabilities) to go from location A (row) to B (column)
#' @param init.structure which state (i.e. location) the initially infected individuals are located.
#'
#' @keywords internal
##

MatrixSanityChecks <- function(structure.matrix, init.structure, none.at.start=NULL) {
  if (!is.matrix(structure.matrix)) stop("structure.matrix should be a matrix.")
  if (ncol(structure.matrix)!=nrow(structure.matrix)) stop("structure.matrix should have the same number of rows and columns.")
  if (!identical(colnames(structure.matrix),rownames(structure.matrix))) stop("structure.matrix rows and columns should have the same names.")
  if (any(rowSums(structure.matrix) != 1)) stop("structure.matrix rows should sum up to 1.")

  if (is.null(none.at.start) && !(init.structure %in% rownames(structure.matrix))) stop("init.structure should be a state present in structure.matrix.")

  if(!is.null(none.at.start) && none.at.start==FALSE && !(init.structure %in% rownames(structure.matrix))) stop("init.structure should be a state present in structure.matrix.")
}

#' @title Checks if the raster is properly formatted
#'
#' @description
#' Checks if the environmental raster is properly formatted given the user's input and nosoi's requirements.
#'
#' @param structure.raster raster object defining the environmental variable.
#' @param init.structure which state (i.e. location) the initially infected individuals are located.
#'
#' @keywords internal
##

RasterSanityChecks <- function(structure.raster, init.structure, none.at.start=NULL) {
  if(!class(structure.raster) == "RasterLayer") stop("structure.raster must be a raster (class RasterLayer).")

  if(is.null(none.at.start)){
    start.env <- raster::extract(structure.raster,cbind(init.structure[1],init.structure[2]))
    if(is.na(start.env)) stop("Your starting position (init.structure) should be on the raster.")
  }

  if(!is.null(none.at.start) && !none.at.start){
    start.env <- raster::extract(structure.raster,cbind(init.structure[1],init.structure[2]))
    if(is.na(start.env)) stop("Your starting position (init.structure) should be on the raster.")
  }

  if(!is.null(none.at.start) && none.at.start){
    start.env <- NA
  }

}

#' @title Progress bar
#'
#' @description
#' Echos the state of the simulation at any given time step provided by the user.
#'
#' @param Host.count.A number of infected hosts of host A.
#' @param Host.count.B number of infected hosts of host B.
#' @param pres.time current time of the simulation
#' @param print.step if print.progress is TRUE, step with which the progress message will be printed.
#' @param length.sim the length (in unit of time) over which the simulation should be run.
#' @param max.infected.A the maximum number of hosts that can be infected in the simulation for host A.
#' @param max.infected.B the maximum number of hosts that can be infected in the simulation for host B.
#' @param type either single/dual host
#'
#' @keywords internal
##

progressMessage <- function(Host.count.A, Host.count.B=NULL, pres.time, print.step, length.sim, max.infected.A, max.infected.B=NULL, type="single") {

  if(type == "single"){
    if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count.A," (",round((Host.count.A/max.infected.A)*100,digits=0),"% of maximum infected hosts).")}
  }

  if(type == "dual"){
    if (pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: (A) ", Host.count.A," (",round((Host.count.A/max.infected.A)*100,digits=0),"% of maximum infected hosts); (B) ", Host.count.B," (",round((Host.count.B/max.infected.B)*100,digits=0),"% of maximum infected hosts).")}
  }
}

#' @title End message
#'
#' @description
#' Message that ends the simulation
#'
#' @param Host.count.A number of infected hosts (host A)
#' @param Host.count.B number of infected hosts (host B)
#' @param pres.time current time of the simulation
#' @param type either single/dual host
#'
#' @keywords internal
##

endMessage <- function(Host.count.A, Host.count.B = NULL, pres.time, type="single") {
  message(endMessageText(Host.count.A, Host.count.B, pres.time, type))
}

endMessageText <- function(Host.count.A, Host.count.B = NULL, pres.time, type="single", done = "done. ") {
  if(type == "single"){
    return(paste0(done, "\nThe simulation has run for ",pres.time," units of time and a total of ",Host.count.A," hosts have been infected."))
  }
  if(type == "dual"){
    return(paste0(done, "\nThe simulation has run for ",pres.time," units of time and a total of ",Host.count.A," (A) and ",Host.count.B, " (B) hosts have been infected."))
  }
}
