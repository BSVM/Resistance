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
#'
#' @examples
#' # Ejemplo con un conjunto de datos de conteos genotipicos para una sola poblacion
#' df1 <- data.frame(
#'  Pop = c("Pop1"),
#'  aa = c(8),
#'  Aa = c(43),
#'  AA = c(25))
#'
#' HWQL1(df1)
#'
#'
#'# Ejemplo de uso con un conjunto de datos de conteos genotípicos para múltiples poblaciones
#'
#'df2 <- data.frame(
#'  Pop = c("Pop 1", "Pop 2", "Pop 3"),
#'  AA = c(4, 0, 1),
#'  Aa = c(33, 39, 22),
#'  aa = c(3, 2, 7))
#'
#'  HWQL1(df2)
#'
#' @export
#'
HWQL1 <- function(genotype_counts) {
  # Input validation
  required_columns <- c("Pop", "AA", "Aa", "aa")
  if (!all(required_columns %in% colnames(genotype_counts))) {
    stop("The data frame 'genotype_counts' must contain the columns: Pop, AA, Aa, aa")
  }
  if (!all(sapply(genotype_counts[, c("AA", "Aa", "aa")], is.numeric))) {
    stop("The columns AA, Aa, and aa must contain numeric values")
  }
  if (any(genotype_counts[, c("AA", "Aa", "aa")] < 0, na.rm = TRUE)) {
    stop("The columns AA, Aa, and aa must contain non-negative values")
  }

  # Initialize the results data frame
  result_df <- data.frame(
    Population = character(0),
    n = numeric(0),
    AA = numeric(0),
    Aa = numeric(0),
    aa = numeric(0),
    Freq_A = numeric(0),
    Freq_a = numeric(0),
    `X2 (df)` = character(0),
    P_Value = numeric(0),
    stringsAsFactors = FALSE
  )

  # Loop to process each row
  for (i in 1:nrow(genotype_counts)) {
    freq_matrix <- genotype_counts[i, , drop = FALSE]

    # Missing value check
    if (any(is.na(freq_matrix[, c("AA", "Aa", "aa")]))) {
      warning(paste("Missing values in population:", genotype_counts$Pop[i]))
      next
    }

    frequencies <- Single_loc_freq(freq_matrix)

    genotypic_freqs <- frequencies[, c("AA", "Aa", "aa")]
    p <- frequencies$Freq_A
    q <- frequencies$Freq_a
    total_individuals <- frequencies$n
    expected_freq_AA <- p^2
    expected_freq_Aa <- 2 * p * q
    expected_freq_aa <- q^2
    expected_counts <- c(
      AA = expected_freq_AA * total_individuals,
      Aa = expected_freq_Aa * total_individuals,
      aa = expected_freq_aa * total_individuals
    )
    observed_counts <- as.numeric(genotype_counts[i, c("AA", "Aa", "aa")])

    # Check that expected values are not zero before calculating chi-squared
    if (any(expected_counts == 0)) {
      warning(paste("Expected count is zero in population:", genotype_counts$Pop[i]))
      next
    }

    chi_squared <- sum((observed_counts - expected_counts)^2 / expected_counts)
    p_value <- pchisq(chi_squared, df = 1, lower.tail = FALSE)

    result_df <- rbind(result_df, data.frame(
      Population = genotype_counts$Pop[i],
      n = total_individuals,
      aa = round(genotypic_freqs$AA, 5),
      Aa = round(genotypic_freqs$Aa, 5),
      AA = round(genotypic_freqs$aa, 5),
      Freq_A = round(p, 5),
      Freq_a = round(q, 5),
      `X2 (df)` = paste0(round(chi_squared, 5), " (1)"),
      P_Value = round(p_value, 5)
    ))
  }

  return(result_df)
}
