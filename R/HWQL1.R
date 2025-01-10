#' Calcular el Test de Hardy-Weinberg
#'
#' Esta función realiza un test de Hardy-Weinberg utilizando los conteos de genotipos observados.
#' Calcula las frecuencias genotípicas y alélicas observadas y esperadas según Hardy-Weinberg,
#' luego evalúa si la población está en equilibrio usando el estadístico de chi-cuadrado y su p-valor.
#'
#' @param genotype_counts Un vector numerico que contiene los conteos de los genotipos observados.
#' Debe tener tres elementos, correspondientes a los genotipos `AA`, `Aa` y `aa` en este orden.
#'
#' @return Un lista con los siguientes elementos:
#' \item{Total_Individuals}{Número total de individuos observados en la muestra.}
#' \item{Genotypic_Frequencies}{Un data frame con las frecuencias genotípicas observadas.}
#' \item{Allelic_Frequencies}{Un data frame con las frecuencias alélicas observadas.}
#' \item{Expected_Frequencies}{Un data frame con las frecuencias esperadas según Hardy-Weinberg para los genotipos `AA`, `Aa` y `aa`, y sus correspondientes conteos esperados.}
#' \item{Observed_Counts}{Un data frame con los conteos observados de los genotipos `AA`, `Aa` y `aa`.}
#' \item{Chi_Squared}{Valor del estadístico chi-cuadrado calculado para la prueba de ajuste Hardy-Weinberg.}
#' \item{P_Value}{P-valor asociado al estadístico chi-cuadrado para determinar si la diferencia entre las frecuencias observadas y esperadas es significativa.}
#'
#' @importFrom stats pchisq
#'
#' @export
#'
#' @examples
#' # Ejemplo con un conjunto de datos de conteos genotípicos
#' mi_tabla <- c(AA = 4, Aa = 33, aa = 3)
#' resultados_HWQ <- calculate_HWQ(mi_tabla)
#' print(resultados_HWQ)
#'
#' # También puedes obtener los conteos observados y las frecuencias esperadas
#' resultados_HWQ$Expected_Frequencies
#' resultados_HWQ$Observed_Counts
#'
calculate_HWQ <- function(genotype_counts) {
  # Calcular frecuencias usando la función anterior
  frequencies <- Single_locus_frequencies(genotype_counts)

  # Extraer frecuencias y el total de individuos
  genotypic_freqs <- frequencies$genotypic_frequencies
  p <- frequencies$allelic_frequencies$Frequency[1]
  q <- frequencies$allelic_frequencies$Frequency[2]
  total_individuals <- frequencies$total_individuals

  # Calcular las frecuencias esperadas según Hardy-Weinberg
  expected_freq_AA <- p^2
  expected_freq_Aa <- 2 * p * q
  expected_freq_aa <- q^2

  # Escalar las frecuencias esperadas a conteos (sin redondeo)
  expected_counts <- c(
    AA = expected_freq_AA * total_individuals,
    Aa = expected_freq_Aa * total_individuals,
    aa = expected_freq_aa * total_individuals
  )

  # Conteos observados
  observed_counts <- genotype_counts

  # Calcular el estadístico de chi-cuadrado
  chi_squared <- sum((observed_counts - expected_counts)^2 / expected_counts)

  # Obtener el p-valor de la prueba chi-cuadrado
  p_value <- pchisq(chi_squared, df = 1, lower.tail = FALSE)

  # Devolver los resultados de la prueba con un formato más limpio
  list(
    Total_Individuals = total_individuals,
    Genotypic_Frequencies = genotypic_freqs,
    Allelic_Frequencies = frequencies$allelic_frequencies,
    Expected_Frequencies = data.frame(
      Genotype = c("AA", "Aa", "aa"),
      Expected_Frequency = c(expected_freq_AA, expected_freq_Aa, expected_freq_aa),
      Expected_Counts = expected_counts
    ),
    Observed_Counts = data.frame(
      Genotype = c("AA", "Aa", "aa"),
      Observed_Count = observed_counts
    ),
    Chi_Squared = chi_squared,
    P_Value = p_value
  )
}

#devtools::load_all()

#devtools::document()
