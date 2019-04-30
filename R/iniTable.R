#' Generates initial table to start the simulation (internal fonction)
#'
#' @description This function creates the initial table for the host.
#'
#' @param init.individuals number of initially infected individuals (i.e. number of lines at time 0).
#' @param init.structure State (or coordinates) of the initially infected individuals.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param ParamHost list of individual based parameters.
#' @param current.environmental.value current value of the environemental variable provided by the raster according to its position in init.structure.
#'
#' @keywords internal

iniTable <- function(init.individuals, init.structure, prefix.host, ParamHost, current.environmental.value=NULL){

  if (init.individuals >= 1){
    list.init <- vector("list", init.individuals)

    for (indiv in 1:init.individuals) {
      list.init[[indiv]] <- newLine(hosts.ID = paste(prefix.host,indiv,sep="-"),
                                    infected.by = paste(NA,indiv,sep="-"),
                                    infected.in = init.structure,
                                    time.is = 0,
                                    ParamHost = ParamHost,
                                    current.environmental.value=current.environmental.value)
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
                                    time.is = 0,
                                    ParamHost = ParamHost,
                                    current.environmental.value=current.environmental.value)
    }
    table.hosts <- data.table::rbindlist(list.init)
    table.hosts <- table.hosts[-1]
  }


  data.table::setkey(table.hosts, "hosts.ID")

  return(table.hosts)
}
