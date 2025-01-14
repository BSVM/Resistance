#' Calculate Allele Frequencies for a Biallelic Marker
#'
#' This function calculates haplotype frequencies for a biallelic genetic marker based on a genotype distribution
#' for multiple populations. The input is a data frame with the genotype counts for each population.
#'
#' @param data A data frame where each row represents a population, the first column is the population name (`Pop`),
#' and the remaining columns are genotype counts (`AABB`, `AABb`, `AAbb`, `AaBB`, `AaBb`, `Aabb`, `aaBB`, `aaBb`, `aabb`).
#'
#' @return A data frame with haplotype frequencies (`AB`, `Ab`, `aB`, `ab`) for each population.
#'
#' @examples
#' # Create a data frame with genotype counts for three populations
#' DicGtps <- data.frame(
#'   Pop = c("Pop 1", "Pop 2", "Pop 3"),
#'   AABB = c(5, 3, 4),
#'   AABb = c(3, 2, 5),
#'   AAbb = c(4, 2, 5),
#'   AaBB = c(4, 6, 3),
#'   AaBb = c(10, 8, 12),
#'   Aabb = c(6, 5, 7),
#'   aaBB = c(3, 4, 2),
#'   aaBb = c(2, 1, 3),
#'   aabb = c(3, 2, 4)
#' )
#'
#' # Calculate haplotype frequencies for the data frame
#' Diloc_Freq_Hd(DicGtps)
#'
#' @export
Diloc_Freq_Hd <- function(data) {
  # Verificar que la columna "Pop" exista
  if (!"Pop" %in% colnames(data)) {
    stop("La columna 'Pop' no existe en el data frame.")
  }

  # Verificar que no haya valores NA en el data frame
  if (any(is.na(data))) {
    stop("El data frame contiene valores NA o faltantes.")
  }

  # Verificar que las columnas necesarias estén presentes
  required_columns <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")
  if (!all(required_columns %in% colnames(data))) {
    stop("El data frame no contiene todas las columnas necesarias.")
  }

  # Verificar que las columnas de los genotipos sean numéricas y no negativas
  if (any(sapply(data[required_columns], function(x) !is.numeric(x) || any(x < 0)))) {
    stop("Las columnas de los genotipos deben ser numéricas y no contener valores negativos.")
  }

  # Inicializar una lista para almacenar los resultados
  resultados <- list()

  # Identificar las columnas de los genotipos
  N_c <- required_columns

  # Iterar sobre cada fila (población)
  for (i in 1:nrow(data)) {
    # Sumar las frecuencias de los genotipos específicos en cada población
    N <- sum(data[i, N_c])

    # Verificar que el número total de individuos no sea cero
    if (N == 0) {
      stop(paste("El número total de individuos en la población", data$Pop[i], "es cero."))
    }

    # Calcular las frecuencias haplotípicas
    f_AB <- (2 * data$AABB[i] + data$AABb[i] + data$AaBB[i] + 0.5 * data$AaBb[i]) / (2 * N)
    f_Ab <- (2 * data$AAbb[i] + data$AABb[i] + data$Aabb[i] + 0.5 * data$AaBb[i]) / (2 * N)
    f_aB <- (2 * data$aaBB[i] + data$AaBB[i] + data$aaBb[i] + 0.5 * data$AaBb[i]) / (2 * N)
    f_ab <- (2 * data$aabb[i] + data$Aabb[i] + data$aaBb[i] + 0.5 * data$AaBb[i]) / (2 * N)

    # Almacenar los resultados en la lista
    resultados[[data$Pop[i]]] <- c(f_AB, f_Ab, f_aB, f_ab)
  }

  # Devolver los resultados como un data frame
  return(data.frame(
    Pop = data$Pop,
    AB = sapply(resultados, function(x) x[1]),
    Ab = sapply(resultados, function(x) x[2]),
    aB = sapply(resultados, function(x) x[3]),
    ab = sapply(resultados, function(x) x[4])
  ))
}
