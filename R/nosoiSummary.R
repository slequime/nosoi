#' Summary of a simulation
#'
#' @description This function provides summary informations about the simulation (number of infected hosts, R0, etc.).
#'
#' @param nosoi.output Output of a nosoi simulation (the data.table).

#' @details All elements are provided in a list.
#'
#' @export nosoiSummary
#' @import stringr
#' @import magrittr
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise

nosoiSummary <- function(nosoi.output){
  nosoiSum = list()
  #number of hosts infected
  nosoiSum[["n.infected"]] = nrow(nosoi.output)

  #number of hosts inactive (have done their full cycle)
  Inactive = nosoi.output[nosoi.output[["active"]] == 0]
  nosoiSum[["n.inactive"]] = nrow(Inactive)

  #estimating R0 (mean number of secondary cases)
  Sec.cases1 = nosoi.output[!str_detect(nosoi.output[["inf.by"]],"NA") & nosoi.output[["inf.by"]] %in% Inactive[["hosts.ID"]]] %>% group_by(inf.by) %>% summarise(Secondary.cases=length(hosts.ID))
  Sec.cases2 = data.frame(inf.by=Inactive[["hosts.ID"]][!Inactive[["hosts.ID"]] %in% Sec.cases1$inf.by],Secondary.cases=0)

  Sec.cases = rbind(Sec.cases,Sec.cases2)

  nosoiSum[["R.node"]] = mean(Sec.cases$Secondary.cases)
  return(nosoiSum)
}
