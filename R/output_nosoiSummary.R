#This file holds all the functions related to parsing & dealing with the output of nosoi (epidemiological side).

#' Summary of a simulation
#'
#' @description This function provides summary informations about the simulation (number of infected hosts, R0, etc.).
#'
#' @param nosoi.output Output of a nosoi simulation (the data.table).

#' @details All elements are provided in a list.
#'
#' @export nosoiSummary
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate

nosoiSummary <- function(nosoi.output){

  #Get R0
  R0 <- getR0(nosoi.output)

  #Get Dynamics
  Dynamics <- getDynamic(nosoi.output)

  #Get Cumulative
  Cumulative <- getCumulative(nosoi.output)

  summary.nosoi = list(R0=R0,
                       dynamics=Dynamics,
                       cumulative=Cumulative)

  return(summary.nosoi)
}

#' @title Number of Active infected units at time t
#'
#' @description
#' For a given time t, this function returns the number of infected active units.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#' @param t time (integer)
#'
#' @return Number of infected units at time t
#'
#' @seealso \code{\link{nosoiSim}}
##
numberInfected <- function(table.nosoi, t) {
  # Attention: need to test that t < t_max
  already_infected <- table.nosoi$inf.time < t  # Strict  inequality: if infected at time t, becomes active at time t + 1
  still_infected <- table.nosoi$out.time > t    # Strict inequality: if out at time t, not active for generation t
  still_infected[is.na(still_infected)] <- TRUE
  return(sum(already_infected & still_infected))
}

#' @title Number of Infected units at time t (BGW)
#'
#' @description
#' For a given time t, this function returns the number of infected active units.
#' The difference with \code{\link{numberInfected}} is that it counts an "out"
#' individual as still there, but with no children.
#' This is for comparision with the BGW process, should be internal only.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#' @param t time (integer)
#'
#' @return Number of infected units at time t
#'
#' @seealso \code{\link{nosoiSim}}
#'
#' @keywords internal
##
numberInfectedBGW <- function(table.nosoi, t) {
  already_infected <- table.nosoi$inf.time <= t  # Non strict inequality: if infected at time t, counts as a member of generation t
  still_infected <- table.nosoi$out.time > t    # Strict inequality: if out at time t, not active for generation t
  still_infected[is.na(still_infected)] <- TRUE
  return(sum(already_infected & still_infected))
}

#' @title Cumulative number of infected units at time t
#'
#' @description
#' For a given time t, this function returns the cumulative number of infected units since the start of the simulation.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#' @param t time (integer)
#'
#' @return Number of infected units at time t
#'
#' @seealso \code{\link{nosoiSim}}
##
cumulativeInfected  <- function(table.nosoi, t) {
  already_infected <- table.nosoi$inf.time < t  # Strict  inequality: if infected at time t, becomes active at time t + 1
  return(sum(already_infected))
}

#' @title Cumulative number of infected units for the full length of the simulation
#'
#' @description
#' Number of infected units at time t for the full length of the simulation
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#'
#'
#' @seealso \code{\link{nosoiSim}}
#' @export getCumulative

getCumulative  <- function(table.nosoi) {

  if (table.nosoi$type == "single"){
    results.cumulative <- data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      temp <- list(t=t, Count = cumulativeInfected(table.nosoi$host.info.A$table.hosts, t),type=table.nosoi$host.info.A$prefix.host)
      results.cumulative <- rbindlist(list(results.cumulative,temp))
    }
  }

  #dual
  if (table.nosoi$type == "dual"){
    results.cumulative <- data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      temp <- list(t=t, Count = cumulativeInfected(table.nosoi$host.info.A$table.hosts, t),type=table.nosoi$host.info.A$prefix.host)
      tempB <- list(t=t, Count = cumulativeInfected(table.nosoi$host.info.B$table.hosts, t),type=table.nosoi$host.info.B$prefix.host)
      results.cumulative <- rbindlist(list(results.cumulative,temp,tempB))
    }
  }
  return(results.cumulative)
}

#' @title Dynamic number of infected units at time t during the full length of the simulation
#'
#' @description
#' Number of infected units at time t for the full length of the simulation
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#'
#' @seealso \code{\link{nosoiSim}}
#' @import magrittr
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @export getDynamic

getDynamic  <- function(table.nosoi) {
  #get rid of some notes linked to the use of dplyr
  time.from <- NULL
  time.to <- NULL
  state <- NULL
  hosts.ID <- NULL

  #get type and pop structure
  #single

  if (table.nosoi$type == "single" & (table.nosoi$host.info.A$popStructure == "none" | table.nosoi$host.info.A$popStructure == "continuous")){
    results.dynamic=data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      temp <- list(t=t, Count = numberInfected(table.nosoi$host.info.A$table.hosts, t),type=table.nosoi$host.info.A$prefix.host)
      results.dynamic <- rbindlist(list(results.dynamic,temp))
    }
  }

  if (table.nosoi$type == "single" & table.nosoi$host.info.A$popStructure == "discrete"){
    results.dynamic=data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      table.state.temp <- subset(table.nosoi$host.info.A$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=table.nosoi$host.info.A$prefix.host,  t=t ) %>% data.table()

      results.dynamic <- rbindlist(list(results.dynamic,table.state.temp))
    }
  }

  #dual
  if (table.nosoi$type == "dual" && (table.nosoi$host.info.A$popStructure == "none" | table.nosoi$host.info.A$popStructure == "continuous") &&
      (table.nosoi$host.info.B$popStructure == "none" || table.nosoi$host.info.B$popStructure == "continuous")){
    results.dynamic=data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      temp <- list(t=t, Count = numberInfected(table.nosoi$host.info.A$table.hosts, t),type=table.nosoi$host.info.A$prefix.host)
      tempB <- list(t=t, Count = numberInfected(table.nosoi$host.info.B$table.hosts, t),type=table.nosoi$host.info.B$prefix.host)
      results.dynamic <- rbindlist(list(results.dynamic,temp,tempB))
    }
  }

  if (table.nosoi$type == "dual" && table.nosoi$host.info.A$popStructure == "discrete" && table.nosoi$host.info.B$popStructure == "discrete"){
    results.dynamic=data.table()
    for(t in 0:(table.nosoi$total.time +1)){
      table.state.tempA <- subset(table.nosoi$host.info.A$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% dplyr::group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=table.nosoi$host.info.A$prefix.host,  t=t ) %>% data.table()

      table.state.tempB <- subset(table.nosoi$host.info.B$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% dplyr::group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=table.nosoi$host.info.B$prefix.host,  t=t ) %>% data.table()

      results.dynamic <- rbindlist(list(results.dynamic,table.state.tempA,table.state.tempB))
    }
  }

  return(results.dynamic)
}

#' @title Get R0 from a nosoi simulation
#'
#' @description
#' Gets an estimate of secondary cases (what R0 usually tries to estimate) and its distribution. The calculation is based on inactive hosts at the end of the simulation.
#'
#' @param table.nosoi result of function \code{\link{nosoiSim}}
#'
#' @return R0 estimate (mean and distribution), as well as the number of inactive hosts at the end of the simulation.
#'
#' @seealso \code{\link{nosoiSim}}
#' @import stringr
#' @import magrittr
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr summarise
#' @export getR0

getR0  <- function(table.nosoi) {
  #To avoids notes (use of dplyr functions)
  suffix <- NULL
  inf.by.y <- NULL
  host.type <- NULL
  hosts.ID <- NULL
  active <- NULL

  if(table.nosoi$type == "single"){
    output.full <- table.nosoi$host.info.A$table.hosts[,c("hosts.ID", "inf.by","active")]
    Inactive <- output.full[output.full[["active"]] == 0]  #get inactive hosts (have done their full cycle)
    n.Inactive <- nrow(Inactive)

    #Secondary case, same host type
    output.small <- output.full[,c("hosts.ID", "inf.by")]
    output.full.merged <- output.full %>% left_join(output.small, by=c("inf.by" = "hosts.ID"), suffix(".x",".y")) %>% as.data.table()

    #estimating R0 (mean number of secondary cases)
    Sec.cases1 <- output.full.merged[!str_detect(output.full.merged[["inf.by"]],"NA") & output.full.merged[["inf.by"]] %in% Inactive[["hosts.ID"]]] %>% group_by(inf.by.y) %>% summarise(Secondary.cases=length(hosts.ID)) %>% data.table()
    Sec.cases2 <- data.table(inf.by.y=as.character(Inactive[["hosts.ID"]][!Inactive[["hosts.ID"]] %in% Sec.cases1$inf.by.y]),Secondary.cases=0) %>% left_join(output.full[,c("hosts.ID")], by=c("inf.by.y"="hosts.ID")) %>% data.table()

    Sec.cases <- rbindlist(list(Sec.cases1,Sec.cases2),use.names=TRUE)
    Sec.cases.A <- Sec.cases$Secondary.cases

    return(list(N.inactive=nrow(Inactive),
                R0.mean=mean(Sec.cases.A),
                R0.dist=Sec.cases.A))
  }

  if(table.nosoi$type == "dual"){
    outputA <- table.nosoi$host.info.A$table.hosts[,c("hosts.ID", "inf.by","active")]
    outputA$host.type <- table.nosoi$host.info.A$prefix.host
    outputB <- table.nosoi$host.info.B$table.hosts[,c("hosts.ID", "inf.by","active")]
    outputB$host.type <- table.nosoi$host.info.B$prefix.host

    output.full = rbindlist(list(outputA,outputB))

    #number of hosts inactive (have done their full cycle)
    Inactive = output.full[output.full[["active"]] == 0]

    N.inactive.A <- nrow(subset(table.nosoi$host.info.A$table.hosts,active==0))
    N.inactive.B <- nrow(subset(table.nosoi$host.info.B$table.hosts,active==0))

    #Secondary case, same host type
    output.small <- output.full[,c("hosts.ID", "inf.by")]
    output.full.merged <- output.full %>% left_join(output.small, by=c("inf.by" = "hosts.ID"), suffix(".x",".y")) %>% as.data.table()

    #estimating R0 (mean number of secondary cases), R0 to the other host
    Sec.cases1 <- output.full.merged[!str_detect(output.full.merged[["inf.by"]],"NA") & output.full.merged[["inf.by"]] %in% Inactive[["hosts.ID"]]] %>% group_by(inf.by.y, host.type) %>% summarise(Secondary.cases=length(hosts.ID)) %>% data.table()
    Sec.cases2 <- data.table(inf.by.y=as.character(Inactive[["hosts.ID"]][!Inactive[["hosts.ID"]] %in% Sec.cases1$inf.by.y]),Secondary.cases=0) %>% left_join(output.full[,c("hosts.ID","host.type")], by=c("inf.by.y"="hosts.ID")) %>% data.table()

    Sec.cases <- rbindlist(list(Sec.cases1,Sec.cases2),use.names=TRUE)

    Sec.cases.A <- subset(Sec.cases, host.type == table.nosoi$host.info.A$prefix.host)$Secondary.cases
    Sec.cases.B <- subset(Sec.cases, host.type == table.nosoi$host.info.B$prefix.host)$Secondary.cases

    return(list(N.inactive.A=N.inactive.A,
                R0.hostA.mean=mean(Sec.cases.A),
                R0.hostA.dist=Sec.cases.A,
                N.inactive.B=N.inactive.B,
                R0.hostB.mean=mean(Sec.cases.B),
                R0.hostB.dist=Sec.cases.B))
  }
}
