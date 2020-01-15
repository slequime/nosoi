#This file holds all the functions related to parsing & dealing with the output of nosoi (epidemiological side).

#' @title Summarizes the epidemiological features of a \code{nosoi} simulation
#'
#' @description This function provides summary information about the simulation (number of infected hosts, R0, etc.) as a list.
#'
#' @param object Output of a nosoi simulation (object of class \code{\link{nosoiSim}}).
#' @param ... further arguments passed to or from other methods.

#' @return All computed data is provided in a list:
#' \describe{
#'    \item{R0}{Provides a sublist with number of inactive hosts at the end of the simulation \code{N.inactive}, mean R0 \code{R0.mean}, and R0 distribution \code{R0.dist}. For more details, see \code{\link{getR0}}.}
#'    \item{dynamics}{\code{\link[data.table:data.table-package]{data.table}} with the count of currently infected (i.e. active) hosts at each time step of the simulation (by state if the simulation was in a discrete structured host population). For more details, see \code{\link{getDynamic}}.}
#'    \item{cumulative}{\code{\link[data.table:data.table-package]{data.table}} with the cumulative count of infected hosts at each time step of the simulation. For more details, see \code{\link{getCumulative}}.}
#'    }
#'
#' @seealso You can directly compute each elements of the list without using the summarise function. See \code{\link{getR0}}, \code{\link{getDynamic}} and \code{\link{getCumulative}}.
#'
#' @examples
#' \donttest{
#'t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#'p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#'p_Exit_fct  <- function(t){return(0.08)}
#'
#'proba <- function(t,p_max,t_incub){
#'  if(t <= t_incub){p=0}
#'  if(t >= t_incub){p=p_max}
#'  return(p)
#'}
#'
#'time_contact <- function(t){round(rnorm(1, 3, 1), 0)}
#'
#'test.nosoi <- nosoiSim(type="single", popStructure="none",
#'                       length=40,
#'                       max.infected=100,
#'                       init.individuals=1,
#'                       nContact=time_contact,
#'                       param.nContact=NA,
#'                       pTrans = proba,
#'                       param.pTrans = list(p_max=p_max_fct,
#'                                           t_incub=t_incub_fct),
#'                       pExit=p_Exit_fct,
#'                       param.pExit=NA)
#'
#'
#' nosoiSummary(test.nosoi)
#'}
#'
#' @export nosoiSummary
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate

nosoiSummary <- function(object){

  #Get R0
  R0 <- getR0(object)

  #Get Dynamics
  Dynamics <- getDynamic(object)

  #Get Cumulative
  Cumulative <- getCumulative(object)

  summary.nosoi = list(R0 = R0,
                       dynamics = Dynamics,
                       cumulative = Cumulative)

  return(summary.nosoi)
}

##
#' @rdname nosoiSummary
#' @export
#' @method summary nosoiSim
##
summary.nosoiSim <- function(object, ...){
  return(nosoiSummary(object))
}

#' @title Number of active infected hosts at time t
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
#'
#' @keywords internal
##
numberInfected <- function(table.nosoi, t) {
  # Attention: need to test that t < t_max
  already_infected <- table.nosoi$inf.time < t  # Strict  inequality: if infected at time t, becomes active at time t + 1
  still_infected <- table.nosoi$out.time > t    # Strict inequality: if out at time t, not active for generation t
  still_infected[is.na(still_infected)] <- TRUE
  return(sum(already_infected & still_infected))
}

#' @title Number of infected hosts at time t (BGW)
#'
#' @description
#' For a given time t, this function returns the number of infected active units.
#' The difference with \code{\link{numberInfected}} is that it counts an "out"
#' individual as still there, but with no children.
#' This is for comparison with the BGW process, should be internal only.
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

#' @title Cumulative number of infected hosts at time t
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
#'
#' @keywords internal
##
cumulativeInfected  <- function(table.nosoi, t) {
  already_infected <- table.nosoi$inf.time < t  # Strict  inequality: if infected at time t, becomes active at time t + 1
  return(sum(already_infected))
}

#' @title Gets the cumulative number of infected hosts for the full length of the simulation
#'
#' @description This function computes from the output of a \code{nosoiSim} simulation the cumulative count of infected hosts at each time step of the simulation. The output is a \code{\link[data.table:data.table-package]{data.table}}.
#'
#' @param nosoi.output Output of a nosoi simulation (object of class \code{\link{nosoiSim}}).
#'
#' @return The output is a \code{\link[data.table:data.table-package]{data.table}} with the following structure:
#' \describe{
#'    \item{t}{Time-step (integer).}
#'    \item{Count}{Cumulative number of infected hosts at given time-step.}
#'    \item{type}{Host-type, identified by its user-defined prefix.}
#'    }
#'
#' @seealso \code{\link{summary.nosoiSim}}
#'
#' @export getCumulative

getCumulative  <- function(nosoi.output) {

  if (nosoi.output$type == "single"){
    results.cumulative <- data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      temp <- list(t=t, Count = cumulativeInfected(nosoi.output$host.info.A$table.hosts, t),type=nosoi.output$host.info.A$prefix.host)
      results.cumulative <- rbindlist(list(results.cumulative,temp))
    }
  }

  #dual
  if (nosoi.output$type == "dual"){
    results.cumulative <- data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      temp <- list(t=t, Count = cumulativeInfected(nosoi.output$host.info.A$table.hosts, t),type=nosoi.output$host.info.A$prefix.host)
      tempB <- list(t=t, Count = cumulativeInfected(nosoi.output$host.info.B$table.hosts, t),type=nosoi.output$host.info.B$prefix.host)
      results.cumulative <- rbindlist(list(results.cumulative,temp,tempB))
    }
  }
  return(results.cumulative)
}

#' @title Gets the current number of infected hosts for the full length of the simulation
#'
#' @description This function computes from the output of a \code{nosoiSim} simulation the dynamic count of infected hosts at each time step (and each state if discrete structure) of the simulation. The output is a \code{\link[data.table:data.table-package]{data.table}}.
#'
#' @param nosoi.output Output of a nosoi simulation (object of class \code{\link{nosoiSim}}).
#'
#' @return The output is a \code{\link[data.table:data.table-package]{data.table}} with the following structure:
#' \describe{
#'    \item{state}{(only when discrete structure) Given state}
#'    \item{Count}{Current number of infected hosts at given time-step.}
#'    \item{type}{Host-type, identified by its user-defined prefix.}
#'    \item{t}{Time-step (integer).}
#'    }
#'
#' @seealso \code{\link{summary.nosoiSim}}
#'
#' @import magrittr
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @export getDynamic

getDynamic  <- function(nosoi.output) {
  #get rid of some notes linked to the use of dplyr
  time.from <- NULL
  time.to <- NULL
  state <- NULL
  hosts.ID <- NULL

  #get type and pop structure
  #single

  if (nosoi.output$type == "single" & (nosoi.output$host.info.A$popStructure == "none" | nosoi.output$host.info.A$popStructure == "continuous")){
    results.dynamic=data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      temp <- list(t=t, Count = numberInfected(nosoi.output$host.info.A$table.hosts, t),type=nosoi.output$host.info.A$prefix.host)
      results.dynamic <- rbindlist(list(results.dynamic,temp))
    }
  }

  if (nosoi.output$type == "single" & nosoi.output$host.info.A$popStructure == "discrete"){
    results.dynamic=data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      table.state.temp <- subset(nosoi.output$host.info.A$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=nosoi.output$host.info.A$prefix.host,  t=t ) %>% data.table()

      results.dynamic <- rbindlist(list(results.dynamic,table.state.temp))
    }
  }

  #dual
  if (nosoi.output$type == "dual" && (nosoi.output$host.info.A$popStructure == "none" | nosoi.output$host.info.A$popStructure == "continuous") &&
      (nosoi.output$host.info.B$popStructure == "none" || nosoi.output$host.info.B$popStructure == "continuous")){
    results.dynamic=data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      temp <- list(t=t, Count = numberInfected(nosoi.output$host.info.A$table.hosts, t),type=nosoi.output$host.info.A$prefix.host)
      tempB <- list(t=t, Count = numberInfected(nosoi.output$host.info.B$table.hosts, t),type=nosoi.output$host.info.B$prefix.host)
      results.dynamic <- rbindlist(list(results.dynamic,temp,tempB))
    }
  }

  if (nosoi.output$type == "dual" && nosoi.output$host.info.A$popStructure == "discrete" && nosoi.output$host.info.B$popStructure == "discrete"){
    results.dynamic=data.table()
    for(t in 0:(nosoi.output$total.time +1)){
      table.state.tempA <- subset(nosoi.output$host.info.A$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% dplyr::group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=nosoi.output$host.info.A$prefix.host,  t=t ) %>% data.table()

      table.state.tempB <- subset(nosoi.output$host.info.B$table.state, (time.from < t & (is.na(time.to)|time.to >= t))) %>% dplyr::group_by(state) %>% dplyr::summarise(Count=length(hosts.ID)) %>% dplyr::mutate(type=nosoi.output$host.info.B$prefix.host,  t=t ) %>% data.table()

      results.dynamic <- rbindlist(list(results.dynamic,table.state.tempA,table.state.tempB))
    }
  }

  return(results.dynamic)
}

#' @title Gets R0 from a \code{nosoi} simulation
#'
#' @description Gets an estimate of secondary cases (what R0 usually tries to estimate) and its distribution from the output of a \code{nosoiSim} simulation. The actual calculation is based on inactive hosts at the end of the simulation to avoid bias introduced by hosts that have not finished their transmission potential.
#'
#' @param nosoi.output Output of a nosoi simulation (object of class \code{\link{nosoiSim}}).
#'
#' @return A list with the following items:
#' \describe{
#'    \item{N.inactive}{Number of inactive hosts at the end of the simulation.}
#'    \item{R0.mean}{Mean R0 based on the distribution (see below).}
#'    \item{R0.dist}{Distribution for each host of the secondary cases it generated (in case of dual-hosts, then the secondary cases of the same host-type).}
#'    }
#'
#' @seealso \code{\link{summary.nosoiSim}}
#'
#' @import stringr
#' @import magrittr
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr summarise
#' @export getR0

getR0  <- function(nosoi.output) {
  #To avoids notes (use of dplyr functions)
  suffix <- NULL
  inf.by.y <- NULL
  host.type <- NULL
  hosts.ID <- NULL
  active <- NULL

  if(nosoi.output$type == "single"){
    output.full <- nosoi.output$host.info.A$table.hosts[,c("hosts.ID", "inf.by","active")]
    Inactive <- output.full[output.full[["active"]] == 0]  #get inactive hosts (have done their full cycle)
    n.Inactive <- nrow(Inactive)
    Sec.cases.A <- NA
    if(n.Inactive > 0) {
    #Secondary case, same host type
    output.small <- output.full[,c("hosts.ID", "inf.by")]
    output.full.merged <- output.full %>% left_join(output.small, by=c("inf.by" = "hosts.ID"), suffix(".x",".y")) %>% as.data.table()

    #estimating R0 (mean number of secondary cases)
    Sec.cases1 <- output.full.merged[!str_detect(output.full.merged[["inf.by"]],"NA") & output.full.merged[["inf.by"]] %in% Inactive[["hosts.ID"]]] %>% group_by(inf.by.y) %>% summarise(Secondary.cases=length(hosts.ID)) %>% data.table()
    Sec.cases2 <- data.table(inf.by.y=as.character(Inactive[["hosts.ID"]][!Inactive[["hosts.ID"]] %in% Sec.cases1$inf.by.y]),Secondary.cases=0) %>% left_join(output.full[,c("hosts.ID")], by=c("inf.by.y"="hosts.ID")) %>% data.table()

    Sec.cases <- rbindlist(list(Sec.cases1,Sec.cases2),use.names=TRUE)
    Sec.cases.A <- Sec.cases$Secondary.cases
    }

    return(list(N.inactive=nrow(Inactive),
                R0.mean=mean(Sec.cases.A),
                R0.dist=Sec.cases.A))
  }

  if(nosoi.output$type == "dual"){
    outputA <- nosoi.output$host.info.A$table.hosts[,c("hosts.ID", "inf.by","active")]
    outputA$host.type <- nosoi.output$host.info.A$prefix.host
    outputB <- nosoi.output$host.info.B$table.hosts[,c("hosts.ID", "inf.by","active")]
    outputB$host.type <- nosoi.output$host.info.B$prefix.host

    output.full = rbindlist(list(outputA,outputB))

    #number of hosts inactive (have done their full cycle)
    Inactive = output.full[output.full[["active"]] == 0]

    N.inactive.A <- nrow(subset(nosoi.output$host.info.A$table.hosts,active==0))
    N.inactive.B <- nrow(subset(nosoi.output$host.info.B$table.hosts,active==0))

    #Secondary case, same host type
    output.small <- output.full[,c("hosts.ID", "inf.by")]
    output.full.merged <- output.full %>% left_join(output.small, by=c("inf.by" = "hosts.ID"), suffix(".x",".y")) %>% as.data.table()

    #estimating R0 (mean number of secondary cases), R0 to the other host
    Sec.cases1 <- output.full.merged[!str_detect(output.full.merged[["inf.by"]],"NA") & output.full.merged[["inf.by"]] %in% Inactive[["hosts.ID"]]] %>% group_by(inf.by.y, host.type) %>% summarise(Secondary.cases=length(hosts.ID)) %>% data.table()
    Sec.cases2 <- data.table(inf.by.y=as.character(Inactive[["hosts.ID"]][!Inactive[["hosts.ID"]] %in% Sec.cases1$inf.by.y]),Secondary.cases=0) %>% left_join(output.full[,c("hosts.ID","host.type")], by=c("inf.by.y"="hosts.ID")) %>% data.table()

    Sec.cases <- rbindlist(list(Sec.cases1,Sec.cases2),use.names=TRUE)

    Sec.cases.A <- subset(Sec.cases, host.type == nosoi.output$host.info.A$prefix.host)$Secondary.cases
    Sec.cases.B <- subset(Sec.cases, host.type == nosoi.output$host.info.B$prefix.host)$Secondary.cases

    return(list(N.inactive.A=N.inactive.A,
                R0.hostA.mean=mean(Sec.cases.A),
                R0.hostA.dist=Sec.cases.A,
                N.inactive.B=N.inactive.B,
                R0.hostB.mean=mean(Sec.cases.B),
                R0.hostB.dist=Sec.cases.B))
  }
}
