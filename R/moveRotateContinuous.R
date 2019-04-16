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
