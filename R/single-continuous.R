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
#' @param structure.raster raster object defining the environmental variable.
#' @param diff.pMove is pMove different between states of the structured population (TRUE/FALSE)
#' @param timeDep.pMove is pMove dependant on the absolute time of the simulation (TRUE/FALSE)
#' @param pMove function that gives the probability of a host moving as a function of time.
#' @param param.pMove parameter names (list of functions) for the pMove.
#' @param diff.moveDist is moveDist dependant on the environmental value (TRUE/FALSE).
#' @param timeDep.moveDist is moveDist dependant on the absolute time of the simulation (TRUE/FALSE).
#' @param moveDist function that gives the distance travelled (based on coordinates); gives the sd value for the brownian motion.
#' @param param.moveDist parameter names (list of functions) for moveDist.
#' @param attracted.by.raster should the hosts be attracted by high values in the environmental raster? (TRUE/FALSE)
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
#' @param progress.bar if TRUE, print message on simulation state to screen.
#' @param print.step progress.bar is TRUE, step with which the progress message will be printed.
#' @param ... other arguments to be passed on to the simulator (see below).
#'
#' @export singleContinuous

singleContinuous <- function(type,
                           structure,
                           length.sim,
                           max.infected,
                           init.individuals,
                           init.structure,
                           structure.raster,
                           diff.pMove=FALSE,
                           timeDep.pMove=FALSE,
                           pMove,
                           param.pMove,
                           diff.moveDist=FALSE,
                           timeDep.moveDist=FALSE,
                           moveDist,
                           param.moveDist,
                           attracted.by.raster=FALSE,
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
  timeContactParsed <- parseFunction(timeContact, param.timeContact, as.character(quote(timeContact)),diff=diff.timeContact, timeDep = timeDep.timeContact, continuous=TRUE)

  #Parsing moveDist
  moveDistParsed <- parseFunction(moveDist, param.moveDist, as.character(quote(moveDist)),diff=diff.moveDist, timeDep = timeDep.moveDist, continuous=TRUE)

  #Parsing pTrans
  pTransParsed <- parseFunction(pTrans, param.pTrans, as.character(quote(pTrans)),diff=diff.pTrans, timeDep = timeDep.pTrans, continuous=TRUE)

  #Parsing pExit
  pExitParsed <- parseFunction(pExit, param.pExit, as.character(quote(pExit)),diff=diff.pExit, timeDep = timeDep.pExit, continuous=TRUE)

  #Discrete states sanity checks -------------------------------------------------------------------------------------------------------------------

  #Extract environmental value at origin:
  RasterSanityChecks(structure.raster,init.structure)
  start.env <- raster::extract(structure.raster,cbind(init.structure[1],init.structure[2]))
  max.raster <- max(structure.raster[], na.rm=T)

  #Parse pMove (same as pExit !!attention if diff)
  pMoveParsed <- parseFunction(pMove, param.pMove, as.character(quote(pMove)),diff=diff.pMove, timeDep = timeDep.pMove, continuous=TRUE)

  #START OF THE SIMULATION --------------------------------------------------------------------------------------------------------

  # Init
  message("Starting the simulation\nInitializing ...", appendLF = FALSE)

  #Creation of initial data ----------------------------------------------------------

  table.hosts <- iniTable(init.individuals, init.structure, prefix.host, param.pExit, param.pMove,param.timeContact, param.pTrans,param.moveDist,current.environmental.value=start.env)

  state.archive <- iniTableState(init.individuals, init.structure, prefix.host,current.environmental.value=start.env)

  Host.count <- init.individuals

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    active.hosts <- table.hosts[["active"]] == 1 #active hosts (boolean vector)
    if (any(active.hosts)){

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
      # state.archive[hosts.ID %in% Move.ID & is.na(time.to), `:=` (time.to = as.numeric(pres.time))]

      state.archive[state.archive[["hosts.ID"]] %in% Move.ID & is.na(state.archive[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

      table.state.temp <- vector("list", length(Move.ID))

      fun <- function(z) {
        moveDistParsed$vect(prestime = pres.time, z[, moveDistParsed$vectArgs, with = FALSE])
      }

      moveDist.values <- table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

      for (i in 1:length(Move.ID)) {

        current.move.pos.x = table.hosts[Move.ID[i],"current.in.x"]
        current.move.pos.y = table.hosts[Move.ID[i],"current.in.y"]
        current.env.value = table.hosts[Move.ID[i],"current.env.value"]
        current.moveDist.value = as.numeric(moveDist.values[i])

        positionFound1 = FALSE
        while (positionFound1 == FALSE)
        {
          counter = 0
          dX = rnorm(1, 0, current.moveDist.value)
          dY = rnorm(1, 0, current.moveDist.value)
          positionFound2 = FALSE
          while (positionFound2 == FALSE)
          {
            angle = (2*base::pi)*runif(1)
            newP = moveRotateContinuous(c(as.numeric(current.move.pos.x),as.numeric(current.move.pos.y)), dX, dY,angle)

            temp.env.value = raster::extract(structure.raster,cbind(newP[1],newP[2]))

            if (!is.na(temp.env.value)){

              if (attracted.by.raster==FALSE) {
                table.hosts[Move.ID[i], `:=` (current.in.x = newP[1])]
                table.hosts[Move.ID[i], `:=` (current.in.y = newP[2])]
                table.hosts[Move.ID[i], `:=` (current.env.value = temp.env.value)]

                table.state.temp[[i]] <- newLineState(Move.ID[i],newP,pres.time,current.environmental.value=temp.env.value)

                positionFound2 = TRUE
                positionFound1 = TRUE
              }

              if (attracted.by.raster==TRUE) {
                counter = counter+1
                v2 = temp.env.value/max.raster
                if (runif(1,0,1) < v2)
                {
                  table.hosts[Move.ID[i], `:=` (current.in.x = newP[1])]
                  table.hosts[Move.ID[i], `:=` (current.in.y = newP[2])]
                  table.hosts[Move.ID[i], `:=` (current.env.value = temp.env.value)]

                  table.state.temp[[i]] <- newLineState(Move.ID[i],newP,pres.time,current.environmental.value=temp.env.value)

                  positionFound2 = TRUE
                  positionFound1 = TRUE
                  if (counter == 30) {
                    positionFound1 = TRUE
                    positionFound2 = TRUE
                  }
                }
              }
            }
          }
        }
      }

      state.archive <- data.table::rbindlist(c(list(state.archive),table.state.temp))
      data.table::setkey(state.archive, "hosts.ID")
    }

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    df.meetTransmit <- table.hosts[active.hosts, c("hosts.ID","current.in.x","current.in.y","current.env.value")]
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

            table.temp[[i]] <- newLine(hosts.ID, as.character(df.meetTransmit[i,]$active.hosts),c(df.meetTransmit[i,]$current.in.x,df.meetTransmit[i,]$current.in.y), pres.time, param.pExit, param.pMove,param.timeContact, param.pTrans,param.moveDist,current.environmental.value=df.meetTransmit[i,]$current.env.value)
            table.state.temp[[i]] <- newLineState(hosts.ID,c(df.meetTransmit[i,]$current.in.x,df.meetTransmit[i,]$current.in.y),pres.time,current.environmental.value=df.meetTransmit[i,]$current.env.value)
          }

          table.hosts <- data.table::rbindlist(c(list(table.hosts),table.temp))
          data.table::setkey(table.hosts,hosts.ID)
          state.archive <- data.table::rbindlist(c(list(state.archive),table.state.temp))
          data.table::setkey(state.archive, "hosts.ID")
        }
      }
    }

    if (progress.bar == TRUE) progressMessage(Host.count, pres.time, print.step, length.sim, max.infected)
    if (Host.count > max.infected) {break}
  }

  endMessage(Host.count, pres.time)

  nosoi.output <- outputWrapper(Host.count, pres.time, table.hosts, state.archive)

  return(nosoi.output)
}
