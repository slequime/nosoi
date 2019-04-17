#' @title Creates a new line to be added to the movement table when hosts moves (internal fonction)
#'
#' @description
#' This function creates a new line for the table,
#' The lines are to be bounded with \code{\link[data.table]{rbindlist}}.
#'
#' @param hosts.ID unique ID for the new host
#' @param state.pres state in which host currently is
#' @param time.is time in the simulation, when the infection takes place
#' @param current.environmental.value current value of environemental variable (from raster) according to coordinates in current.in.
#'
#' @return a list with the new line to add.
#'
#' @keywords internal

newLineState <- function(hosts.ID,state.pres,time.is,current.environmental.value=NA) {

  if (length(state.pres) == 1){
    return(list(hosts.ID = hosts.ID,
                state = state.pres,
                time.from = time.is,
                time.to = NA_real_
    )
    )
  }

  if (length(state.pres) == 2){
  return(list(hosts.ID = hosts.ID,
           state.x = state.pres[1],
           state.y = state.pres[2],
           current.env.value = current.environmental.value,
           time.from = time.is,
           time.to = NA_real_
  )
  )
}
}
