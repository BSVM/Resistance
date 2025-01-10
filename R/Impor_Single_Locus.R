#' Importa datos de un solo locus
#'
#' Esta funcion verifica la consistencia de una matriz de genotipos y cuenta el numero de individuos con cada genotipo (aa, Aa, AA) en un conjunto de datos.
#'
#' @param df Un data.frame llamado "Matriz de Genotipos" que representa la presencia de los genotipos (aa, Aa, AA) para un conjunto de individuos identificados por su ID.
#'
#' @return Un vector numerico con el conteo de individuos para cada genotipo en el orden: AA, Aa, aa.
#'
#' @examples
#'
#'
#' genotipos <- data.frame(
#'   ID = 1:12,
#'   aa = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
#'   Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
#'   AA = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Impor_Single_Locus(genotipos)
#'
#' @export
Impor_Single_Locus <- function(df) {
  # Verificar si hay mas de una "X" en una misma fila
  for (i in 1:nrow(df)) {
    if (sum(df[i, c("aa", "Aa", "AA")] == "X") > 1) {
      stop(paste("Error: Mas de una 'X' en la fila", i))
    }
  }

  # Contar el numero total de "X"
  total_X <- sum(df == "X")

  # Verificar si el conteo de "X" es mayor o menor al numero de filas
  if (total_X > nrow(df) || total_X < nrow(df)) {
    stop("Error: El conteo de 'X' es mayor o menor al numero de filas")
  }

  # Contar el numero de individuos con cada genotipo
  aa <- sum(df$aa == "X")
  Aa <- sum(df$Aa == "X")
  AA <- sum(df$AA == "X")

  # Crear la tabla de frecuencias
  mi_tabla <- as.numeric(c(AA, Aa, aa))

  return(mi_tabla)
}
