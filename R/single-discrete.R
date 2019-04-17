#' Single-host with structured host population
#'
#' @description This function runs a single-host transmission chain simulation, with a structured host population (such as spatial features).
#' The simulation stops either at the end of given time (specified by length.sim) or when the number of hosts infected threshold (max.infected)
#' is passed.
#'
#' @param type single host or dual host
#' @param structure is there any discrete structure in the host population to be accounted for (such as spatial features)
#' @param length.sim specifies the length (in unit of time) over which the simulation should be run.
#' @param max.infected specifies the maximum number of hosts that can be infected in the simulation.
#' @param init.individuals number of initially infected individuals.
#' @param init.structure which state (i.e. location) the initially infected individuals are located.
#' @param structure.matrix transition matrix (probabilities) to go from location A (row) to B (column)
#' @param diff.pMove is pMove different between states of the structured population (TRUE/FALSE)
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.timeContact is timeContact different between states of the structured population (TRUE/FALSE)
#' @param timeContact function that gives the number of potential transmission events per unit of time.
#' @param param.timeContact parameter names (list of functions) for timeContact.
#' @param diff.pTrans is pTrans different between states of the structured population (TRUE/FALSE)
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pTrans.
#' @param diff.pExit is pExit different between states of the structured population (TRUE/FALSE)
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param progress.bar if TRUE, displays a progress bar (current time/length.sim).
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @export singleDiscrete

singleDiscrete <- function(type,
                           structure,
                           length.sim,
                           max.infected,
                           init.individuals,
                           init.structure,
                           structure.matrix,
                           diff.pMove=FALSE,
                           timeDep.pMove=FALSE,
                           pMove,
                           param.pMove,
                           diff.timeContact=FALSE,
                           timeDep.timeContact=FALSE,
                           timeContact,
                           param.timeContact,
                           diff.pTrans=FALSE,
                           timeDep.pTrans=FALSE,
                           pTrans,
                           param.pTrans,
                           diff.pExit=FALSE,
                           timeDep.pExit=FALSE,
                           pExit,
                           param.pExit,
                           prefix.host="H",
                           progress.bar=TRUE,
                           print.step=10,
                           ...){

  #Sanity checks---------------------------------------------------------------------------------------------------------------------------
  #This section checks if the arguments of the function are in a correct format for the function to run properly

  if (is.na(length.sim) | length.sim <= 1) stop("You must specify a length (in time units) for your simulation (bigger than 1).")
  if (is.na(max.infected) | max.infected <= 1) stop("You must specify a maximum number of infected hosts (bigger than 1).")
  if (is.na(init.individuals) | init.individuals < 1 | !init.individuals%%1==0) stop("The transmission chain should be started by 1 or more (integer) individuals")

  #Parsing timeContact
  if (diff.timeContact == FALSE & timeDep.timeContact == FALSE) {
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)))
  }

  if (diff.timeContact == TRUE & timeDep.timeContact == FALSE) {
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),'current.in'))==FALSE) stop("timeContact should have 'current.in' as a variable. diff.timeContact == TRUE.")
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("timeContact should have a realisation for each possible state. diff.timeContact is TRUE.")
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=TRUE)
  }

  if (diff.timeContact == FALSE & timeDep.timeContact == TRUE) {
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),'prestime'))==FALSE) stop("timeContact should have 'prestime' as a variable. timeDep.timeContact == TRUE.")
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=FALSE,timeDep=TRUE)
  }

  if (diff.timeContact == TRUE & timeDep.timeContact == TRUE) {
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),'current.in'))==FALSE) stop("timeContact should have 'current.in' as a variable. diff.timeContact == TRUE.")
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("timeContact should have a realisation for each possible state. diff.timeContact is TRUE.")
    if (any(str_detect(paste0(as.character(body(timeContact)),collapse=" "),'prestime'))==FALSE) stop("timeContact should have 'prestime' as a variable. timeDep.timeContact == TRUE.")
    timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=TRUE,timeDep=TRUE)
  }

  #Parsing pTrans
  if (diff.pTrans == FALSE & timeDep.pTrans == FALSE) {
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)))
  }

  if (diff.pTrans == TRUE & timeDep.pTrans == FALSE) {
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),'current.in'))==FALSE) stop("pTrans should have 'current.in' as a variable. diff.pTrans == TRUE.")
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pTrans should have a realisation for each possible state. diff.pTrans is TRUE.")
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=TRUE)
  }

  if (diff.pTrans == FALSE & timeDep.pTrans == TRUE) {
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),'prestime'))==FALSE) stop("pTrans should have 'prestime' as a variable. timeDep.pTrans == TRUE.")
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=FALSE,timeDep=TRUE)
  }

  if (diff.pTrans == TRUE & timeDep.pTrans == TRUE) {
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),'current.in'))==FALSE) stop("pTrans should have 'current.in' as a variable. diff.pTrans == TRUE.")
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pTrans should have a realisation for each possible state. diff.pTrans is TRUE.")
    if (any(str_detect(paste0(as.character(body(pTrans)),collapse=" "),'prestime'))==FALSE) stop("pTrans should have 'prestime' as a variable. timeDep.pTrans == TRUE.")
    pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=TRUE,timeDep=TRUE)
  }

  #Parsing pExit

  if (diff.pExit == FALSE & timeDep.pExit == FALSE) {
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)))
  }

  if (diff.pExit == TRUE & timeDep.pExit == FALSE) {
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),'current.in'))==FALSE) stop("pExit should have 'current.in' as a variable. diff.pExit == TRUE.")
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pExit should have a realisation for each possible state. diff.pExit is TRUE.")
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=TRUE)
  }

  if (diff.pExit == FALSE & timeDep.pExit == TRUE) {
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),'prestime'))==FALSE) stop("pExit should have 'prestime' as a variable. timeDep.pExit == TRUE.")
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=FALSE,timeDep=TRUE)
  }

  if (diff.pExit == TRUE & timeDep.pExit == TRUE) {
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),'current.in'))==FALSE) stop("pExit should have 'current.in' as a variable. diff.pExit == TRUE.")
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pExit should have a realisation for each possible state. diff.pExit is TRUE.")
    if (any(str_detect(paste0(as.character(body(pExit)),collapse=" "),'prestime'))==FALSE) stop("pExit should have 'prestime' as a variable. timeDep.pExit == TRUE.")
    pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=TRUE,timeDep=TRUE)
  }

  #Discrete states sanity checks -------------------------------------------------------------------------------------------------------------------

  if (!is.matrix(structure.matrix)) stop("structure.matrix should be a matrix.")
  if (ncol(structure.matrix)!=nrow(structure.matrix)) stop("structure.matrix should have the same number of rows and columns.")
  if (!identical(colnames(structure.matrix),rownames(structure.matrix))) stop("structure.matrix rows and columns should have the same names.")
  if (any(rowSums(structure.matrix) != 1)) stop("structure.matrix rows should sum up to 1.")
  if (!init.structure %in% rownames(structure.matrix)) stop("init.structure should be a state present in structure.matrix.")

  melted.structure.matrix <- reshape2::melt(structure.matrix, varnames = c("from","to"),value.name="prob", as.is = TRUE) #melting the matrix go get from -> to in one line with probability

  #Parse pMove (same as pExit !!attention if diff)

  if (diff.pMove == FALSE & timeDep.pMove == FALSE) {
    pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)))
  }

  if (diff.pMove == TRUE & timeDep.pMove == FALSE) {
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),'current.in'))==FALSE) stop("pMove should have 'current.in' as a variable. diff.pMove == TRUE.")
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pMove should have a realisation for each possible state. diff.pMove is TRUE.")
    pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=TRUE)
  }

  if (diff.pMove == FALSE & timeDep.pMove == TRUE) {
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),'prestime'))==FALSE) stop("pMove should have 'prestime' as a variable. timeDep.pMove == TRUE.")
    pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=FALSE,timeDep=TRUE)
  }

  if (diff.pMove == TRUE & timeDep.pMove == TRUE) {
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),'current.in'))==FALSE) stop("pMove should have 'current.in' as a variable. diff.pMove == TRUE.")
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),paste0('current.in == "',colnames(structure.matrix),'"'))==FALSE)) stop("pMove should have a realisation for each possible state. diff.pMove is TRUE.")
    if (any(str_detect(paste0(as.character(body(pMove)),collapse=" "),'prestime'))==FALSE) stop("pMove should have 'prestime' as a variable. timeDep.pMove == TRUE.")
    pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=TRUE,timeDep=TRUE)
  }

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation")
  message("Initializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  table.hosts <- iniTable(init.individuals, init.structure, prefix.host, param.pExit, param.pMove,param.timeContact, param.pTrans,param.moveDist=NA)

  state.archive <- iniTableState(init.individuals, init.structure, prefix.host)

  Host.count <- init.individuals

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    active.hosts <- table.hosts[["active"]] == 1 #active hosts (boolean vector)

    if (any(active.hosts)) {

      fun <- function(z) {
        pExitParsed$vect(prestime = pres.time, z[, pExitParsed$vectArgs, with = FALSE])
      }
      p.exit.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

      exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
    }

    exiting.full <- active.hosts
    exiting.full[exiting.full] <- exiting

    table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                    active = 0)]

    exiting.ID <- table.hosts[active.hosts][as.vector(exiting), "hosts.ID"]$hosts.ID

    state.archive[state.archive[["hosts.ID"]] %in% exiting.ID & is.na(state.archive[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

    active.hosts[active.hosts] <- !exiting # Update active hosts

    if (!any(active.hosts)) {break} #if no more active hosts, then end simulation

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    fun <- function(z) {
      pMoveParsed$vect(prestime = pres.time, z[, pMoveParsed$vectArgs, with = FALSE])
    }
    p.move.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]
    moving <- drawBernouilli(p.move.values) #Draws K bernouillis with various probability (see function for more detail)
    # }

    moving.full <- active.hosts
    moving.full[moving.full] <- moving

    #step 1.2 if moving, where are they going?

    Move.ID <- table.hosts[moving.full,][["hosts.ID"]]

    if (length(Move.ID) > 0){
      #Updating state archive for moving individuals:

      state.archive[state.archive[["hosts.ID"]] %in% Move.ID & is.na(state.archive[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

      table.state.temp <- vector("list", length(Move.ID))

      for (i in 1:length(Move.ID)) {

        current.move.pos <- melted.structure.matrix[which(melted.structure.matrix$from==as.character(table.hosts[Move.ID[i],"current.in"])),]

        going.to <- sample(current.move.pos$to, 1, replace = FALSE, prob = current.move.pos$prob)
        table.state.temp[[i]] <- newLineState(Move.ID[i],going.to,pres.time)
        table.hosts[Move.ID[i], `:=` (current.in = going.to)]

      }

      state.archive <- data.table::rbindlist(c(list(state.archive),table.state.temp))
      data.table::setkey(state.archive, "hosts.ID")
    }

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    df.meetTransmit <- table.hosts[active.hosts, c("hosts.ID","current.in")]
    df.meetTransmit[, active.hosts:=hosts.ID]

    fun <- function(z) {
      timeContactParsed$vect(prestime = pres.time, z[, timeContactParsed$vectArgs, with = FALSE])
    }
    timeContact.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

    df.meetTransmit$number.contacts <- timeContact.values

    haveContact <- df.meetTransmit[["number.contacts"]] > 0
    df.meetTransmit <- df.meetTransmit[haveContact]
    active.hosts[active.hosts] <- haveContact # Update active hosts

    if (nrow(df.meetTransmit) > 0) {

      fun <- function(z) {
        pTransParsed$vect(prestime = pres.time, z[, pTransParsed$vectArgs, with = FALSE])
      }

      df.meetTransmit[, "Ptransmit"] <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"] #adds transmission probability to events
      df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Ptransmit"]] > 0] #discards event with probability 0

      if (nrow(df.meetTransmit) > 0) {

        df.meetTransmit <- df.meetTransmit[rep(seq(1, nrow(df.meetTransmit)), df.meetTransmit$number.contacts)]

        df.meetTransmit[,"Trans"] <- drawBernouilli(df.meetTransmit[["Ptransmit"]]) #Draws K bernouillis with various probability (see function for more detail)

        df.meetTransmit <- df.meetTransmit[df.meetTransmit[["Trans"]]] #Discards events with no realisation

        if (nrow(df.meetTransmit) >0) {
          table.temp <- vector("list", nrow(df.meetTransmit))
          table.state.temp <- vector("list", nrow(df.meetTransmit))
          for (i in 1:nrow(df.meetTransmit)) {

            Host.count <- Host.count+1
            hosts.ID <- as.character(paste(prefix.host,Host.count,sep="-"))

            table.temp[[i]] <- newLine(hosts.ID, as.character(df.meetTransmit[i,]$active.hosts),as.character(df.meetTransmit[i,]$current.in), pres.time, param.pExit, param.pMove,param.timeContact, param.pTrans,param.moveDist=NA)
            table.state.temp[[i]] <- newLineState(hosts.ID,as.character(df.meetTransmit[i,]$current.in),pres.time)
          }

          table.hosts <- data.table::rbindlist(c(list(table.hosts),table.temp))
          data.table::setkey(table.hosts,hosts.ID)
          state.archive <- data.table::rbindlist(c(list(state.archive),table.state.temp))
          data.table::setkey(state.archive, "hosts.ID")
        }
      }
    }

    if (progress.bar == TRUE & pres.time%%print.step == 0) {message("Time: ", pres.time ," (",round((pres.time/length.sim)*100,digits=0),"% of maximum length). Hosts count: ", Host.count," (",round((Host.count/max.infected)*100,digits=0),"% of maximum infected hosts).")}
    if (Host.count > max.infected) {break}
  }
  message(" done.")
  message("The simulation has run for ",pres.time," units of time and a total of ",Host.count," hosts have been infected.")

  nosoi.output <- list()

  nosoi.output[["total.time"]] <- pres.time
  nosoi.output[["N.infected"]] <- Host.count
  nosoi.output[["table.hosts"]] <- table.hosts
  nosoi.output[["table.state"]] <- state.archive

  return(nosoi.output)
}
