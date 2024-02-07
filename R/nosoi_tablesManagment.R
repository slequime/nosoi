#This file holds all the function needed to create and update the tables used in nosoi

#' @title Creates a new line to be added to the table when new host is infected (internal function)
#'
#' @description
#' This function creates a new line for the table.
#' The lines are to be bounded with \code{\link[data.table]{rbindlist}}.
#'
#' @param hosts.ID unique ID for the new host
#' @param infected.by unique ID of host that transmits to the new one
#' @param infected.in state in which the host was infected
#' @param time.is time in the simulation, when the infection takes place
#' @param ParamHost list of individual based parameters.
#' @param current.environmental.value current environmental value
#' @param current.cell.number.raster unique number of the raster cell where the host is
#' @param current.count.A current count of host A
#' @param current.count.B current count of host B
#'
#' @return a list with the new line to add.
#'
#' @keywords internal

newLine <- function(hosts.ID,
                    infected.by,
                    infected.in,
                    time.is,
                    ParamHost,
                    current.environmental.value = NULL,
                    current.cell.number.raster = NULL,
                    current.count.A = integer(0),
                    current.count.B = integer(0)) {

  if (length(infected.in) == 1) {
    if (is.na(infected.in)) infected.in <- NULL
    return(c(hosts.ID = hosts.ID,
             inf.by = infected.by,
             inf.in = infected.in,
             current.in = infected.in,
             host.count.A = as.integer(current.count.A),
             host.count.B = as.integer(current.count.B),
             inf.time = time.is,
             out.time = NA_integer_,
             active = TRUE,
             as.list(sapply(ParamHost, function(x) x(1)))
    )
    )
  }

  if (length(infected.in) == 2) {

    return(c(hosts.ID = hosts.ID,
             inf.by = infected.by,
             inf.in.x = infected.in[1],
             inf.in.y = infected.in[2],
             current.in.x = infected.in[1],
             current.in.y = infected.in[2],
             current.env.value = current.environmental.value,
             current.cell.raster = current.cell.number.raster,
             host.count.A = as.integer(current.count.A),
             host.count.B = as.integer(current.count.B),
             inf.time = time.is,
             out.time = NA_integer_,
             active = TRUE,
             as.list(sapply(ParamHost, function(x) x(1)))
    )
    )
  }
}

#' Generates initial table to start the simulation (internal function)
#'
#' @description This function creates the initial table for the host.
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param init.structure State (or coordinates) of the initially infected individuals.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param ParamHost list of individual based parameters.
#' @param current.environmental.value current value of the environmental  variable provided by the raster according to its position in init.structure.
#' @param current.cell.number.raster unique number of the raster cell where the host is
#' @param current.count.A current count of host A
#' @param current.count.B current count of host B
#'
#' @keywords internal

iniTable <- function(init.individuals, init.structure, prefix.host, ParamHost,
                     current.environmental.value = NULL, current.cell.number.raster = NULL,
                     current.count.A = integer(0), current.count.B = integer(0)){

  if (init.individuals >= 1){
    list.init <- vector("list", init.individuals)

    for (indiv in 1:init.individuals) {
      list.init[[indiv]] <- newLine(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                    infected.by = paste(NA,indiv,sep="-"),
                                    infected.in = init.structure,
                                    time.is = 0L,
                                    ParamHost = ParamHost,
                                    current.environmental.value = current.environmental.value,
                                    current.cell.number.raster = current.cell.number.raster,
                                    current.count.A = current.count.A,
                                    current.count.B = current.count.B)
    }
    table.hosts <- data.table::rbindlist(list.init)
  }

  if (init.individuals == 0){
    fake.init.individuals <- 1
    list.init <- vector("list", fake.init.individuals)

    for (indiv in 1:fake.init.individuals) {
      list.init[[indiv]] <- newLine(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                    infected.by = paste(NA,indiv,sep="-"),
                                    infected.in = init.structure,
                                    time.is = 0L,
                                    ParamHost = ParamHost,
                                    current.environmental.value = current.environmental.value,
                                    current.cell.number.raster = current.cell.number.raster,
                                    current.count.A = current.count.A,
                                    current.count.B = current.count.B)
    }
    table.hosts <- data.table::rbindlist(list.init)
    table.hosts <- table.hosts[-1]
  }


  data.table::setkey(table.hosts, "hosts.ID")

  return(table.hosts)
}

#' @title Creates a new line to be added to the movement table when hosts moves (internal function)
#'
#' @description
#' This function creates a new line for the table,
#' The lines are to be bounded with \code{\link[data.table]{rbindlist}}.
#'
#' @param hosts.ID unique ID for the new host
#' @param state.pres state in which host currently is
#' @param time.is time in the simulation, when the infection takes place
#' @param current.environmental.value current value of environmental  variable (from raster) according to coordinates in current.in.
#' @param current.cell.number.raster unique number of the raster cell where the host is
#' @return a list with the new line to add.
#'
#' @keywords internal

newLineState <- function(hosts.ID, state.pres, time.is, current.environmental.value = NA, current.cell.number.raster = NA) {

  if (length(state.pres) == 1){
    return(list(hosts.ID = hosts.ID,
                state = state.pres,
                time.from = time.is,
                time.to = NA_integer_
    )
    )
  }

  if (length(state.pres) == 2){
    return(list(hosts.ID = hosts.ID,
                state.x = state.pres[1],
                state.y = state.pres[2],
                current.env.value = current.environmental.value,
                current.cell.raster = current.cell.number.raster,
                time.from = time.is,
                time.to = NA_integer_
    )
    )
  }
}

#' Generates initial movement table to start the simulation (internal function)
#'
#' @description This function creates the initial table for the host, with 5+number of parameters of the transmission probability function parameters, and init.individuals row(s).
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param init.structure State of the initially infected individuals.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param current.environmental.value current value of the environmental  variable provided by the raster according to its position in init.structure.
#' @param current.cell.number.raster unique number of the raster cell where the host is
#'
#' @keywords internal

iniTableState <- function(init.individuals, init.structure, prefix.host, current.environmental.value = NULL, current.cell.number.raster = NA){

  if (init.individuals >= 1){
    list.init <- vector("list", init.individuals)

    for (indiv in 1:init.individuals) {
      list.init[[indiv]] <- newLineState(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                         state.pres = init.structure,
                                         time.is = 0L,
                                         current.environmental.value = current.environmental.value,
                                         current.cell.number.raster = current.cell.number.raster)
    }

    table.mov <- data.table::rbindlist(list.init)

  }

  if (init.individuals == 0){
    fake.init.individuals <- 1
    list.init <- vector("list", fake.init.individuals)

    for (indiv in 1:fake.init.individuals) {
      list.init[[indiv]] <- newLineState(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                         state.pres = init.structure,
                                         time.is = 0L,
                                         current.environmental.value=current.environmental.value,
                                         current.cell.number.raster=current.cell.number.raster)
    }
    table.mov <- data.table::rbindlist(list.init)
    table.mov <- table.mov[-1]
  }
  data.table::setkey(table.mov, "hosts.ID")

  return(table.mov)
}

#' @title Update table state with exiting individuals
#'
#' @description
#' Return a vector of exiting individuals.
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @param exiting a vector (TRUE/FALSE) of exiting individuals
#' @param pres.time current time
#'
#' @return \code{nosoiSim} object
#'
#' @keywords internal
##

updateTableState <- function(res, exiting, pres.time) {
  #To avoids notes (use of dplyr functions)
  time.to <- NULL
  hosts.ID <- NULL

  exiting.ID <- res$table.hosts[exiting]$hosts.ID
  exiting.state <- res$table.state[,is.na(time.to) & hosts.ID %in% exiting.ID]

  res$table.state[exiting.state, `:=` (time.to = pres.time)]

  return(res)
}
