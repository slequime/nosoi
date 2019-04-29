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

  res <- nosoiSimConstructor(N.infected = init.individuals,
                             pres.time = 1,
                             table.hosts = iniTable(init.individuals, init.structure, prefix.host, param.pExit, param.pMove, param.timeContact, param.pTrans, param.moveDist, current.environmental.value = start.env),
                             table.state = iniTableState(init.individuals, init.structure, prefix.host, current.environmental.value = start.env),
                             type = "singleContinuous")

  # Running the simulation ----------------------------------------
  message(" running ...")

  for (pres.time in 1:length.sim) {

    #Step 0: Active hosts ----------------------------------------------------------
    active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)
    if (any(active.hosts)){

      fun <- function(z) {
        pExitParsed$vect(prestime = pres.time, z[, pExitParsed$vectArgs, with = FALSE])
      }
      p.exit.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][["V1"]]

      exiting <- drawBernouilli(p.exit.values) #Draws K bernouillis with various probability (see function for more detail)
    }

    exiting.full <- active.hosts
    exiting.full[exiting.full] <- exiting

    res$table.hosts[exiting.full, `:=` (out.time = as.numeric(pres.time),
                                        active = 0)]

    exiting.ID <- res$table.hosts[active.hosts][as.vector(exiting), "hosts.ID"]$hosts.ID

    res$table.state[res$table.state[["hosts.ID"]] %in% exiting.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

    active.hosts[active.hosts] <- !exiting # Update active hosts

    if (!any(active.hosts)) {break} #if no more active hosts, then end simulation

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
      # res$table.state[hosts.ID %in% Move.ID & is.na(time.to), `:=` (time.to = as.numeric(pres.time))]

      res$table.state[res$table.state[["hosts.ID"]] %in% Move.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

      table.state.temp <- vector("list", length(Move.ID))

      fun <- function(z) {
        moveDistParsed$vect(prestime = pres.time, z[, moveDistParsed$vectArgs, with = FALSE])
      }

      moveDist.values <- res$table.hosts[active.hosts, fun(.SD), by="hosts.ID"][, "V1"]

      for (i in 1:length(Move.ID)) {

        current.move.pos.x = res$table.hosts[Move.ID[i],"current.in.x"]
        current.move.pos.y = res$table.hosts[Move.ID[i],"current.in.y"]
        current.env.value = res$table.hosts[Move.ID[i],"current.env.value"]
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
                res$table.hosts[Move.ID[i], `:=` (current.in.x = newP[1])]
                res$table.hosts[Move.ID[i], `:=` (current.in.y = newP[2])]
                res$table.hosts[Move.ID[i], `:=` (current.env.value = temp.env.value)]

                table.state.temp[[i]] <- newLineState(Move.ID[i],newP,pres.time,current.environmental.value=temp.env.value)

                positionFound2 = TRUE
                positionFound1 = TRUE
              }

              if (attracted.by.raster==TRUE) {
                counter = counter+1
                v2 = temp.env.value/max.raster
                if (runif(1,0,1) < v2)
                {
                  res$table.hosts[Move.ID[i], `:=` (current.in.x = newP[1])]
                  res$table.hosts[Move.ID[i], `:=` (current.in.y = newP[2])]
                  res$table.hosts[Move.ID[i], `:=` (current.env.value = temp.env.value)]

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

      res$table.state <- data.table::rbindlist(c(list(res$table.state),table.state.temp))
      data.table::setkey(res$table.state, "hosts.ID")
    }

    #Step 2: Hosts Meet & Transmist ----------------------------------------------------

    res <- meetTransmit(res,
                        pres.time,
                        active.hosts,
                        positions = c("current.in.x", "current.in.y", "current.env.value"),
                        timeContactParsed, pTransParsed,
                        prefix.host, param.pExit, param.pMove, param.timeContact, param.pTrans,
                        param.moveDist)

    if (progress.bar == TRUE) progressMessage(res$N.infected, pres.time, print.step, length.sim, max.infected)
    if (res$N.infected > max.infected) {break}
  }

  endMessage(res$N.infected, pres.time)

  res$pres.time <- pres.time

  return(res)
}
