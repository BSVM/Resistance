#' Calcular frecuencias genotípicas y alélicas
#'
#' Esta función calcula las frecuencias genotípicas y alélicas a partir de una matriz que contiene los conteos de genotipos (AA, Aa, aa) por población.
#'
#' @param freq_matrix Un data frame con las siguientes columnas:
#' \describe{
#'   \item{Pop}{Un identificador para cada población (e.g., "Pop1", "Pop2").}
#'   \item{aa}{El conteo de individuos con el genotipo aa.}
#'   \item{Aa}{El conteo de individuos con el genotipo Aa.}
#'   \item{AA}{El conteo de individuos con el genotipo AA.}
#' }
#'
#' @return Un data frame con las siguientes columnas:
#' \describe{
#'   \item{Population}{El identificador de la población.}
#'   \item{n}{El número total de individuos considerados en los cálculos.}
#'   \item{aa}{La frecuencia genotípica de individuos con el genotipo aa.}
#'   \item{Aa}{La frecuencia genotípica de individuos con el genotipo Aa.}
#'   \item{AA}{La frecuencia genotípica de individuos con el genotipo AA.}
#'   \item{Freq_A}{La frecuencia alélica del alelo A.}
#'   \item{Freq_a}{La frecuencia alélica del alelo a.}
#' }
#'
#' @examples
#' # Ejemplo: Matriz de frecuencias genotípicas para una población
#' freq_matrix <- data.frame(
#'   Pop = c("Pop1"),
#'   aa = c(2),
#'   Aa = c(9),
#'   AA = c(1))
#'
#' # Calcular frecuencias
#' Single_loc_freq(freq_matrix)
#'
#'# Ejemplo: Matriz de frecuencias genotípicas para multiples poblaciones
#'
#'
#'   df2 <- Import_single_loc(Loc1Pops3)
#'
#'   Single_loc_freq(df2)
#' @export
Single_loc_freq <- function(freq_matrix) {
  # Crear un dataframe vacío para almacenar los resultados
  result_df <- data.frame(
    Pop = character(0),
    n = integer(0),
    aa = numeric(0),
    Aa = numeric(0),
    AA = numeric(0),
    Freq_A = numeric(0),
    Freq_a = numeric(0),
    stringsAsFactors = FALSE
  )

  # Iterar sobre cada fila de la matriz de frecuencias por población
  for (i in 1:nrow(freq_matrix)) {
    # Extraer los conteos de genotipos para cada población
    genotype_counts <- as.numeric(freq_matrix[i, c("aa", "Aa", "AA")])

    # Calcular el número total de individuos
    total_individuals <- sum(genotype_counts)

    # Calcular las frecuencias genotípicas
    freq_AA <- genotype_counts[3] / total_individuals
    freq_Aa <- genotype_counts[2] / total_individuals
    freq_aa <- genotype_counts[1] / total_individuals

    # Calcular las frecuencias alélicas
    p <- (genotype_counts[3] + 0.5 * genotype_counts[2]) / total_individuals  # Frecuencia del alelo A
    q <- 1 - p  # Frecuencia del alelo a

    # Almacenar los resultados en el dataframe
    result_df <- rbind(result_df, data.frame(
      Pop = freq_matrix$Pop[i],
      n = total_individuals,
      aa = freq_aa,
      Aa = freq_Aa,
      AA = freq_AA,
      Freq_A = p,
      Freq_a = q
    ))
  }

  return(result_df)
}

# devtools::load_all()

# devtools::document()
