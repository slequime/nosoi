#' @title Make Move function
#'
#' @description
#' Makes the move.
#'
#' @param res an object of class \code{nosoiSimOne}.
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
                      sdMoveParsed = NA, structure.raster = NA,
                      attracted.by.raster = NA, max.raster = NA) {

  moveFunction <-  switch(res$popStructure,
                          discrete = moveFunction.discrete,
                          continuous = moveFunction.continuous)

  return(moveFunction(res, pres.time, moving.full,
                      structure.matrix = structure.matrix,
                      sdMoveParsed = sdMoveParsed,
                      structure.raster = structure.raster,
                      attracted.by.raster = attracted.by.raster,
                      max.raster = max.raster))
}

#' @title Make Discrete Move function
#'
#' @description
#' Makes the discrete move.
#'
#' @param res an object of class \code{nosoiSimOne}.
#' @param pres.time current time
#' @inheritParams singleDiscrete
#'
#' @return The modified object res
#'
#' @keywords internal
##
moveFunction.discrete <- function(res, pres.time, moving.full, structure.matrix, ...) {

  Move.ID <- res$table.hosts[moving.full,][["hosts.ID"]]

  #melting the matrix go get from -> to in one line with probability
  melted.structure.matrix <- reshape2::melt(structure.matrix,
                                            varnames = c("from", "to"),
                                            value.name = "prob",
                                            as.is = TRUE)

  if (length(Move.ID) > 0){
    #Updating state archive for moving individuals:

    res$table.state[res$table.state[["hosts.ID"]] %in% Move.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = pres.time)]

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
#' @param res an object of class \code{nosoiSimOne}.
#' @param pres.time current time
#' @inheritParams singleContinuous
#'
#' @return The modified object res
#'
#' @keywords internal
##
moveFunction.continuous <- function(res, pres.time, moving.full,
                                    sdMoveParsed, structure.raster,
                                    attracted.by.raster, max.raster, ...) {

  Move.ID <- res$table.hosts[moving.full,][["hosts.ID"]]
  Move.index <- which(moving.full)

  if (length(Move.ID) > 0){
    #Updating state archive for moving individuals:
    # res$table.state[hosts.ID %in% Move.ID & is.na(time.to), `:=` (time.to = as.numeric(pres.time))]
    active.hosts <- res$table.hosts[["active"]] #active hosts (boolean vector)

    res$table.state[res$table.state[["hosts.ID"]] %in% Move.ID & is.na(res$table.state[["time.to"]]), `:=` (time.to = pres.time)]

    table.state.temp <- vector("list", length(Move.ID))

    sdMove.values <- applyFunctionToHosts(res, pres.time, sdMoveParsed, active.hosts)

    for (i in 1:length(Move.ID)) {

      current.move.pos <- c(res$table.hosts[[Move.index[i], "current.in.x"]], res$table.hosts[[Move.index[i], "current.in.y"]])
      # current.env.value = res$table.hosts[Move.index[i],"current.env.value"]
      current.sdMove.value = sdMove.values[[i]]

      positionFound1 <- FALSE
      counter1 <- 0
      while (positionFound1 == FALSE && counter1 < 30) {
        counter1 <- counter1 + 1

        dX = rnorm(1, 0, current.sdMove.value)
        dY = rnorm(1, 0, current.sdMove.value)

        counter2 <- 0
        positionFound2 = FALSE
        while (positionFound2 == FALSE && counter2 <= 30) { # If counter2 > 30, try another move
          counter2 = counter2 + 1

          angle = (2*base::pi)*runif(1)
          newP = moveRotateContinuous(current.move.pos, dX, dY, angle)

          temp.cell.number = raster::cellFromXY(structure.raster, cbind(newP[1],newP[2]))
          temp.env.value = raster::extract(structure.raster,temp.cell.number)

          if (!is.na(temp.env.value)){

            if (!attracted.by.raster) {
              set(res$table.hosts, Move.index[i],
                  c("current.in.x", "current.in.y", "current.env.value", "current.cell.raster"),
                  list(newP[1], newP[2], temp.env.value, temp.cell.number))

              table.state.temp[[i]] <- newLineState(Move.ID[i],newP,pres.time,current.environmental.value=temp.env.value,
                                                    current.cell.number.raster=temp.cell.number)

              positionFound2 = TRUE
              positionFound1 = TRUE
            }

            if (attracted.by.raster) {
              v2 = temp.env.value/max.raster
              if ((runif(1,0,1) < v2) || (counter2 >= 30)) { # If counter2 == 30, just accept anyway.
                set(res$table.hosts, Move.index[i],
                    c("current.in.x", "current.in.y", "current.env.value", "current.cell.raster"),
                    list(newP[1], newP[2], temp.env.value, temp.cell.number))

                table.state.temp[[i]] <- newLineState(Move.ID[i],newP,pres.time,current.environmental.value=temp.env.value,
                                                      current.cell.number.raster=temp.cell.number)

                positionFound2 = TRUE
                positionFound1 = TRUE
              }
            }
          }
        }
      }
      if (counter1 == 30) { # Just impossible to find a correct move, stay where you are
        table.state.temp[[i]] <- newLineState(Move.ID[i],
                                              current.move.pos, pres.time,
                                              current.environmental.value = res$table.hosts[Move.index[i],"current.env.value"],
                                              current.cell.number.raster = res$table.hosts[Move.index[i],"current.cell.number.raster"])
      }
    }

    res$table.state <- data.table::rbindlist(c(list(res$table.state),table.state.temp))
    data.table::setkey(res$table.state, "hosts.ID")
  }

  return(res)
}

## Rotates movement in 2D space

#' @title  Rotates movement in 2D space
#'
#' @description
#'  Rotates movement in 2D space
#'
#' @param pt1 current position (x,y)
#' @param dX movement in X dim
#' @param dY movement in Y space
#' @param angle angle

#' @keywords internal
##

moveRotateContinuous = function(pt1, dX, dY, angle)
{
  s = sin(angle); c = cos(angle)
  x = dX; y = dY
  x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
  x_new = x_new+pt1[1]; y_new = y_new+pt1[2]
  return(c(x_new,y_new))
}
