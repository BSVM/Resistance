#' calcula minimo y maximo de un vector aleatorio
#'
#' @param n tamano de la muestra.
#' @param media media de la poblacion.
#' @param desvi desviacion de la poblacion.
#'
#' @return un vector con minimo y maximo.
#'
#' @examples
#' min_max(n=10, media=50, desvi=2)
#' min_max(n=100, desvi=3, media=180)
#'
#' @importFrom stats rnorm
#'
#' @export
min_max <- function(n, media, desvi) {
  x <- rnorm(n=n, mean=media, sd=desvi)
  res <- c("minimo"=min(x), "max"=max(x))
  return(res)
}

#
#devtools::load_all()

#devtools::document()

# devtools::install()


