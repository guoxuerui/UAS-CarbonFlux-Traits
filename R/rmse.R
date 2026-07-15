#' Root Mean Squared Error
#' 
#' @param x numerical vector
#' @param y numerical vector
#' @param na.rm remove NA?



rmse <- function(x, y, na.rm = TRUE) {
  sqrt(mean((x - y) ^ 2, na.rm = na.rm))
}