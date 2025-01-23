#' Realiza un Test de Equilibrio de Hardy-Weinberg para múltiples poblaciones con marcadores dihíbridos (dos loci)
#'
#' Esta función calcula las frecuencias genotípicas y haplotípicas observadas y esperadas
#' según el equilibrio de Hardy-Weinberg para múltiples poblaciones con marcadores dihíbridos
#' (dos loci). Utiliza el método de iteración de Expectation-Maximization (EM) para estimar
#' las frecuencias haplotípicas y realiza una prueba de chi-cuadrado para evaluar el equilibrio
#' de Hardy-Weinberg. **Esta función está diseñada específicamente para marcadores dihíbridos**,
#' es decir, dos loci con dos alelos cada uno (A/a y B/b), lo que genera genotipos como AABB,
#' AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb y aabb.
#'
#' @param genotype_counts Un data frame que contiene los conteos de genotipos observados
#'   para múltiples poblaciones. Debe incluir las siguientes columnas:
#'   \itemize{
#'     \item \code{Pop}: Identificador de la población.
#'     \item \code{AABB}, \code{AABb}, \code{AAbb}, \code{AaBB}, \code{AaBb}, \code{Aabb},
#'       \code{aaBB}, \code{aaBb}, \code{aabb}: Conteos de los genotipos para los dos loci.
#'   }
#' @param gl Grados de libertad para la prueba de chi-cuadrado. Por defecto es 1.
#'
#' @return Un data frame con los siguientes resultados para cada población:
#'   \itemize{
#'     \item \code{Population}: Identificador de la población.
#'     \item \code{n}: Número total de individuos en la población.
#'     \item \code{AB}, \code{Ab}, \code{aB}, \code{ab}: Frecuencias haplotípicas estimadas
#'       para los haplotipos AB, Ab, aB y ab.
#'     \item \code{X2 (df)}: Valor del estadístico chi-cuadrado y grados de libertad.
#'     \item \code{P_Value}: P-valor asociado al estadístico chi-cuadrado.
#'   }
#'
#' @details
#' La función utiliza la función \code{Diloc_Freq_Hd} para estimar las frecuencias
#' haplotípicas mediante el método de iteración EM. Luego, calcula las frecuencias
#' genotípicas esperadas bajo el equilibrio de Hardy-Weinberg y realiza una prueba de
#' chi-cuadrado para comparar las frecuencias observadas y esperadas.
#'
#' Si una población tiene menos de dos genotipos con conteos mayores que cero, se omite
#' y se muestra una advertencia.
#'
#' @importFrom stats pchisq
#'
#' @examples
#' # Cargar datos de ejemplo
#' data <- Import_di_loc(Loc2Pops3)
#'
#' # Realizar el test de Hardy-Weinberg con 3 grados de libertad
#' resultados <- HWQL2(data$DicGtps, gl = 3)
#'
#' # Mostrar los resultados
#' print(resultados)
#'
#' @export
HWQL2 <- function(genotype_counts, gl = 1) {
  # Input validation
  required_columns <- c("Pop", "AABB", "AABb", "AAbb","AaBB", "AaBb", "Aabb", "aaBB","aaBb", "aabb")
  if (!all(required_columns %in% colnames(genotype_counts))) {
    stop("The data frame 'genotype_counts' must contain the columns: Pop, AA, Aa, aa")
  }
  if (!all(sapply(genotype_counts[, c("AABB", "AABb", "AAbb","AaBB", "AaBb", "Aabb", "aaBB","aaBb", "aabb")], is.numeric))) {
    stop("The columns AA, Aa, and aa must contain numeric values")
  }
  if (any(genotype_counts[, c("AABB", "AABb", "AAbb","AaBB", "AaBb", "Aabb", "aaBB","aaBb", "aabb")] < 0, na.rm = TRUE)) {
    stop("The columns AA, Aa, and aa must contain non-negative values")
  }

  # Initialize the results data frame
  result_df <- data.frame(
    Population = character(0),
    n = numeric(0),
    AB = numeric(0),
    Ab = numeric(0),
    aB = numeric(0),
    ab = numeric(0),
    `X2 (df)` = character(0),
    P_Value = numeric(0),
    stringsAsFactors = FALSE
  )

  # Loop to process each row
  for (i in 1:nrow(genotype_counts)) {
    freq_matrix <- genotype_counts[i, , drop = FALSE]

    # Missing value check
    if (any(is.na(freq_matrix[, c("AABB", "AABb", "AAbb","AaBB", "AaBb", "Aabb", "aaBB","aaBb", "aabb")]))) {
      warning(paste("Missing values in population:", genotype_counts$Pop[i]))
      next
    }

    # Call Diloc with the current row's data
    df <- Diloc_Freq_Hd(freq_matrix, max_iterations = 10, tol = 1e-6, dat = 1)

    # Check if Diloc returned the expected results
    if (is.null(df$Frecuencias_Haplotipos) || is.null(df$Frecuencias_Haplotipos$AB) ||
        is.null(df$Frecuencias_Haplotipos$Ab) || is.null(df$Frecuencias_Haplotipos$aB) ||
        is.null(df$Frecuencias_Haplotipos$ab) || is.null(df$Frecuencias_Haplotipos$N)) {
      warning(paste("Diloc did not return expected results for population:", genotype_counts$Pop[i]))
      next
    }

    # Extract haplotype frequencies
    AB <- df$Frecuencias_Haplotipos$AB
    Ab <- df$Frecuencias_Haplotipos$Ab
    aB <- df$Frecuencias_Haplotipos$aB
    ab <- df$Frecuencias_Haplotipos$ab
    total_individuals <- df$Frecuencias_Haplotipos$N

    # Calculate expected frequencies
    expected_freq_AABB <- if (freq_matrix$AABB > 0) AB^2 else 0
    expected_freq_AABb <- if (freq_matrix$AABb > 0) 2 * AB * Ab else 0
    expected_freq_AAbb <- if (freq_matrix$AAbb > 0) Ab^2 else 0
    expected_freq_AaBB <- if (freq_matrix$AaBB > 0) 2 * AB * aB else 0
    expected_freq_AaBb <- if (freq_matrix$AaBb > 0) 2 * (AB * ab + Ab * aB) else 0
    expected_freq_Aabb <- if (freq_matrix$Aabb > 0) 2 * Ab * ab else 0
    expected_freq_aaBB <- if (freq_matrix$aaBB > 0) aB^2 else 0
    expected_freq_aaBb <- if (freq_matrix$aaBb > 0) 2 * aB * ab else 0
    expected_freq_aabb <- if (freq_matrix$aabb > 0) ab^2 else 0

    # Calculate expected counts
    expected_counts <- c(
      AABB =  expected_freq_AABB * total_individuals,
      AABb =  expected_freq_AABb * total_individuals,
      AAbb =  expected_freq_AAbb * total_individuals,
      AaBB =  expected_freq_AaBB * total_individuals,
      AaBb =  expected_freq_AaBb * total_individuals,
      Aabb =  expected_freq_Aabb * total_individuals,
      aaBB =  expected_freq_aaBB * total_individuals,
      aaBb =  expected_freq_aaBb * total_individuals,
      aabb =  expected_freq_aabb * total_individuals
    )

    observed_counts <- as.numeric(freq_matrix[, c("AABB", "AABb", "AAbb","AaBB", "AaBb", "Aabb", "aaBB","aaBb", "aabb")])

    # Filter out genotypes with zero expected counts
    valid_genotypes <- expected_counts > 0
    expected_counts <- expected_counts[valid_genotypes]
    observed_counts <- observed_counts[valid_genotypes]

    # Check that there are valid genotypes to calculate chi-squared
    if (length(expected_counts) == 0) {
      warning(paste("No valid genotypes to calculate chi-squared in population:", genotype_counts$Pop[i]))
      next
    }

    chi_squared <- sum((observed_counts - expected_counts)^2 / expected_counts)
    p_value <- pchisq(chi_squared, df = gl, lower.tail = FALSE)  # Usar el valor de gl proporcionado

    # Add results to the data frame
    result_df <- rbind(result_df, data.frame(
      Population = genotype_counts$Pop[i],
      n = total_individuals,
      AB = round(AB, 5),
      Ab = round(Ab, 5),
      aB = round(aB, 5),
      ab = round(ab, 5),
      `X2 (df)` = paste0(round(chi_squared, 5), " (", gl, ")"),  # Mostrar los grados de libertad
      P_Value = round(p_value, 5),
      stringsAsFactors = FALSE
    ))
  }

  return(result_df)
}

