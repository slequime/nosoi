#' Generates initial movement table to start the simulation (internal fonction)
#'
#' @description This function creates the initial table for the host, with 5+number of parameters of the transmission probability function paramters, and init.individuals row(s).
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param init.structure State of the initially infected individuals.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param current.environmental.value current value of the environemental variable provided by the raster according to its position in init.structure.
#'
#' @keywords internal

iniTableState <- function(init.individuals, init.structure, prefix.host, current.environmental.value=NULL){

  list.init <- vector("list", init.individuals)

  for (indiv in 1:init.individuals) {
      list.init[[indiv]] <- newLineState(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                    state.pres = init.structure,
                                    time.is = 0,
                                    current.environmental.value=current.environmental.value)
  }

   table.mov <- data.table::rbindlist(list.init)
   data.table::setkey(table.mov, "hosts.ID")

  return(table.mov)
}
