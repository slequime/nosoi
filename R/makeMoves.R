#' @title Make Move function
#'
#' @description
#' Makes the move.
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @inheritParams singleDiscrete
#' @inheritParams singleContinuous
#'
#' @return The modified object res
#'
#' @keywords internal
##
makeMoves <- function(res, pres.time, moving.full,
                      structure.matrix = NA,
                      moveDistParsed = NA, structure.raster = NA,
                      attracted.by.raster = NA, max.raster = NA) {

  moveFunction <-  switch(res$type,
                          singleDiscrete = moveFunction.singleDiscrete,
                          singleContinuous = moveFunction.singleContinuous)

  Move.ID <- res$table.hosts[moving.full,][["hosts.ID"]]

  return(moveFunction(res, pres.time, Move.ID,
                      structure.matrix = structure.matrix,
                      moveDistParsed = moveDistParsed,
                      structure.raster = structure.raster,
                      attracted.by.raster = attracted.by.raster,
                      max.raster = max.raster))
}

#' @title Make Discrete Move function
#'
#' @description
#' Makes the discrete move.
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @inheritParams singleDiscrete
#'
#' @return The modified object res
#'
#' @keywords internal
##
moveFunction.singleDiscrete <- function(res, pres.time, Move.ID, structure.matrix, ...) {

  #melting the matrix go get from -> to in one line with probability
  melted.structure.matrix <- reshape2::melt(structure.matrix,
                                            varnames = c("from", "to"),
                                            value.name = "prob",
                                            as.is = TRUE)

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

  return(res)
}

#' @title Make Continuous Move function
#'
#' @description
#' Makes the continuous move.
#'
#' @param res an object of class \code{nosoiSim}.
#' @param pres.time current time
#' @inheritParams singleContinuous
#'
#' @return The modified object res
#'
#' @keywords internal
##
moveFunction.singleContinuous <- function(res, pres.time, moving.full,
                                          moveDistParsed, structure.raster,
                                          attracted.by.raster, max.raster, ...) {

  Move.ID <- res$table.hosts[moving.full,][["hosts.ID"]]

  if (length(Move.ID) > 0){
    #Updating state archive for moving individuals:
    # res$table.state[hosts.ID %in% Move.ID & is.na(time.to), `:=` (time.to = as.numeric(pres.time))]
    active.hosts <- res$table.hosts[["active"]] == 1 #active hosts (boolean vector)

    res$table.state[res$table.state[["hosts.ID"]] %in% Move.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = as.numeric(pres.time))]

    table.state.temp <- vector("list", length(Move.ID))

    moveDist.values <- applyFunctionToHosts(res, pres.time, moveDistParsed, active.hosts)

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

  return(res)
}
