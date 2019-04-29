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
#' @param timeDep.pMove is pMove dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.timeContact is timeContact different between states of the structured population (TRUE/FALSE)
#' @param timeDep.timeContact is timeContact dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param timeContact function that gives the number of potential transmission events per unit of time.
#' @param param.timeContact parameter names (list of functions) for timeContact.
#' @param diff.pTrans is pTrans different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pTrans is pTrans dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pTrans function that gives the probability of transmit a pathogen as a function of time since infection.
#' @param param.pTrans parameter names (list of functions) for the pTrans.
#' @param diff.pExit is pExit different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pExit is pExit dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pExit function that gives the probability to exit the simulation for an infected host (either moving out, dying, etc.).
#' @param param.pExit parameter names (list of functions) for the pExit.
#' @param prefix.host character(s) to be used as a prefix for the hosts identification number.
#' @param progress.bar if TRUE, displays a progress bar (current time/length.sim).
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
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

  CoreSanityChecks(length.sim, max.infected, init.individuals)

  #Parsing timeContact
  timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=diff.timeContact, timeDep = timeDep.timeContact, stateNames=colnames(structure.matrix))

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=diff.pTrans, timeDep = timeDep.pTrans, stateNames=colnames(structure.matrix))

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=diff.pExit, timeDep = timeDep.pExit, stateNames=colnames(structure.matrix))

  #Discrete states sanity checks -------------------------------------------------------------------------------------------------------------------

  MatrixSanityChecks(structure.matrix,init.structure)

  melted.structure.matrix <- reshape2::melt(structure.matrix, varnames = c("from","to"),value.name="prob", as.is = TRUE) #melting the matrix go get from -> to in one line with probability

  #Parse pMove (same as pExit !!attention if diff)
  pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=diff.pMove, timeDep = timeDep.pMove, stateNames=colnames(structure.matrix))

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  res <- nosoiSimConstructor(N.infected = init.individuals,
                             total.time = 1,
                             table.hosts = iniTable(init.individuals, init.structure, prefix.host, param.pExit, param.pMove, param.timeContact, param.pTrans, param.moveDist = NA),
                             table.state = iniTableState(init.individuals, init.structure, prefix.host),
                             type = "singleDiscrete")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    exiting.full <- getExiting(res, pres.time, pExitParsed)

    res$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                        active = 0)]

    active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

    if (!any(active.hosts)) {break}

    #Step 1: Moving ----------------------------------------------------

    #step 1.1 which hosts are moving

    fun <- function(z) {
      pMoveParsed$vect(prestime = pres.time, z[, pMoveParsed$vectArgs, with = FALSE])
    }
    p.move.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][["V1"]]
    moving <- drawBernouilli(p.move.values) #Draws K bernouillis with various probability (see function for more detail)
    # }

    moving.full <- active.hosts
    moving.full[moving.full] <- moving

    #step 1.2 if moving, where are they going?

    Move.ID <- res$table.hosts[moving.full,][["hosts.ID"]]

    if (length(Move.ID) > 0){
      #Updating state archive for moving individuals:

      res$table.state[res$table.state[["hosts.ID"]] %in% Move.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

      table.state.temp <- vector("list", length(Move.ID))

      for (i in 1:length(Move.ID)) {

        current.move.pos <- melted.structure.matrix[which(melted.structure.matrix$from==as.character(res$table.hosts[Move.ID[i],"current.in"])),]

        going.to <- sample(current.move.pos$to, 1, replace = FALSE, prob = current.move.pos$prob)
        table.state.temp[[i]] <- newLineState(Move.ID[i],going.to,pres.time)
        res$table.hosts[Move.ID[i], `:=` (current.in = going.to)]

      }

      res$table.state <- data.table::rbindlist(c(list(res$table.state),table.state.temp))
      data.table::setkey(res$table.state, "hosts.ID")
    }

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    res <- meetTransmit(res,
                        pres.time,
                        active.hosts,
                        positions = c("current.in"),
                        timeContactParsed, pTransParsed,
                        prefix.host, param.pExit, param.pMove, param.timeContact, param.pTrans,
                        param.moveDist = NA)

    if (progress.bar == TRUE) progressMessage(res$N.infected, pres.time, print.step, length.sim, max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(res$N.infected, pres.time)

  res$total.time <- pres.time

  return(res)
}
