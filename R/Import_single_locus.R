#' Import Single Locus Data Mutation KDR
#'
#' This function checks the consistency of a genotype matrix and counts the number of individuals with each genotype (aa, Aa, AA) in a dataset. It is specifically designed for working with KDR mutation genotypes that confer resistance to pyrethroids.
#'
#' @param df A data.frame called "Genotype Matrix" that represents the presence of genotypes (aa, Aa, AA) for a set of individuals identified by their ID.
#'
#' @return A numeric vector with the count of individuals for each genotype in the order: AA, Aa, aa.
#'
#' @details
#' The function performs several checks to ensure data integrity:
#' - Validates that the input is a `data.frame`.
#' - Ensures the presence of the required columns: `ID`, `aa`, `Aa`, and `AA`.
#' - Checks for the absence of `NA` values in the required columns.
#' - Ensures that each genotype column only contains the values: `X` (indicating the individual has the genotype) and `0` (indicating absence of the genotype).
#' - Ensures that each individual (row) has exactly one `X`, meaning an individual can only present one genotype.
#' - Validates that the total count of `X` matches the total number of individuals (rows) in the data frame.
#'
#'
#'
#' @examples
#'
#' # Example 1: Provides a matrix indicating the genotype for 12 individuals for 1 population
#'
#' df <- data.frame(
#'   ID = 1:12,
#'   Pop = rep("Pop1", 12),
#'   AA = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
#'   Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
#'   aa = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#'
#' Import_single_loc(df)
#'
#'
#' # Example 2: Provides a matrix indicating the genotype for 106 individuals in 3 populations
#'
#' Import_single_loc(Loc1Pops3)
#'
#' @export
Import_single_loc <- function(df) {
  # Validar el tipo de dato
  if (!is.data.frame(df)) {
    stop("Error: The input must be a data.frame")
  }

  # Verificar el orden de las columnas
  correct_order <- c("ID", "Pop", "AA", "Aa", "aa")
  if (!all(colnames(df) == correct_order)) {
    stop("Error: The columns must be in the following order: ID, Pop, AA, Aa, aa")
  }

  # Verificar si hay NAs en el dataframe
  if (any(is.na(df))) {
    cols_with_na <- colnames(df)[colSums(is.na(df)) > 0]
    stop(paste("Error: The data.frame contains NA values in the following columns:", paste(cols_with_na, collapse=", ")))
  }

  # Verificar caracteres inválidos
  invalid_chars <- df[, 3:5] != "X" & df[, 3:5] != "0"
  if (any(invalid_chars)) {
    cols_with_invalid <- colnames(df)[3:5][apply(invalid_chars, 2, any)]
    stop(paste("Error: The data.frame contains invalid characters in the following columns:", paste(cols_with_invalid, collapse=", "), ". Only 'X' and '0' are allowed"))
  }

  # Verificar si hay más de una "X" en cada fila y especificar la fila y columnas con error
  rows_with_errors <- list()
  for (i in 1:nrow(df)) {
    count_X <- sum(df[i, 3:5] == "X") # Cuenta cuántas "X" hay en las columnas AA, Aa y aa de la fila i
    if (count_X != 1) { # Si no hay exactamente una "X" en la fila
      rows_with_errors[[length(rows_with_errors) + 1]] <- list(row = i, count_X = count_X)
    }
  }

  if (length(rows_with_errors) > 0) {
    error_messages <- sapply(rows_with_errors, function(x) {
      paste("Row", x$row, "has", x$count_X, "'X'")
    })
    stop(paste("Error: The following rows do not have exactly one 'X':\n", paste(error_messages, collapse = "\n")))
  }

  # Contar el número total de "X" en cada columna
  total_X_AA <- sum(df$AA == "X")
  total_X_Aa <- sum(df$Aa == "X")
  total_X_aa <- sum(df$aa == "X")
  total_X <- total_X_AA + total_X_Aa + total_X_aa

  # Verificar si el conteo de 'X' es igual al número de filas
  if (total_X != nrow(df)) {
    # Salida del conteo total de X y número de filas para depuración
    print(paste("Total X count in AA:", total_X_AA))
    print(paste("Total X count in Aa:", total_X_Aa))
    print(paste("Total X count in aa:", total_X_aa))
    print(paste("Total X count:", total_X))
    print(paste("Number of rows:", nrow(df)))

    stop("Error: The count of 'X' must be equal to the number of rows")
  }

  # Crear un data.frame vacío para almacenar las frecuencias por población
  freq_matrix <- data.frame(
    Pop = character(length(unique(df$Pop))),
    AA = numeric(length(unique(df$Pop))),
    Aa = numeric(length(unique(df$Pop))),
    aa = numeric(length(unique(df$Pop))),
    stringsAsFactors = FALSE
  )

  # Obtener poblaciones únicas
  populations <- unique(df$Pop)

  for (i in seq_along(populations)) {
    pop <- populations[i]
    df_pop <- df[df$Pop == pop, ]

    # Contar el número de individuos con cada genotipo para la población
    AA<- sum(df_pop$AA == "X")
    Aa <- sum(df_pop$Aa == "X")
    aa <- sum(df_pop$aa == "X")

    # Almacenar los resultados en el data.frame
    freq_matrix[i, ] <- list(pop, AA, Aa, aa)
  }

  return(freq_matrix)
}
