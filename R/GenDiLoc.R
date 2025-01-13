#' Calcular frecuencia genotípica de un marcador dilocus
#'
#' Esta función calcula las frecuencias genotípicas para un marcador genético dilocus
#' en una población, basado en los genotipos observados en los individuos.
#'
#' @param Data Un data frame que contiene los datos geneticos. Debe incluir al menos
#' las columnas `Pop` (identificador de la población) y las columnas de genotipos:
#' `AA`, `Aa`, `aa`, `BB`, `Bb`, `bb`, que indican la presencia de cada combinación alelica.
#'
#' @return Una matriz en la que cada fila corresponde a una población y cada columna
#' representa un genotipo posible (`AABB`, `AABb`, etc.). Las celdas contienen las frecuencias
#' absolutas (conteos) de cada genotipo en cada población.
#'
#' @export
#'
#' @examples
#' # Not Run sin utilizar la funcion Import_di_loc
#' df <- data.frame(
#'   ID = 1:12,
#'   Pop = rep("Pop1", 12),
#'   AA = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
#'   Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
#'   aa = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'   BB = c("X", 0, 0, 0, 0, 0, "X", 0, 0, 0, 0, 0),
#'   Bb = c(0, "X", "X", "X", 0, "X", 0, 0, "X", "X", "X", "X"),
#'   bb = c(0, 0, 0, 0, "X", 0, 0, "X", 0, 0, 0, 0))
#'
#' # Not Run
#' GenDiLoc(df)
#'
GenDiLoc <- function(Data) {
  # Check for NA values
  if (any(is.na(Data))) {
    stop("Error: The data contains NA values.")
  }

  # Check that all values are "0" or "X" except for the "Pop" and "ID" columns
  valid_values <- c("0", "X")
  for (col in names(Data)[!(names(Data) %in% c("Pop", "ID"))]) {  # Exclude the "Pop" and "ID" columns
    if (!all(Data[[col]] %in% valid_values)) {
      stop("Error: The data contains values other than '0' and 'X'.")
    }
  }

  # Check that the first column is "Pop"
  if (names(Data)[2] != "Pop") {
    stop("Error: The second column must be 'Pop'.")
  }

  # Check that all required columns are present
  required_columns <- c("Pop", "AA", "Aa", "aa", "BB", "Bb", "bb")
  missing_columns <- setdiff(required_columns, names(Data))
  if (length(missing_columns) > 0) {
    stop(paste("Error: The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }

  # Create a new column 'Genotype' based on the conditions
  Data$Genotype <- with(Data, ifelse(AA == "X" & BB == "X", "AABB",
                                     ifelse(AA == "X" & Bb == "X", "AABb",
                                            ifelse(AA == "X" & bb == "X", "AAbb",
                                                   ifelse(Aa == "X" & BB == "X", "AaBB",
                                                          ifelse(Aa == "X" & Bb == "X", "AaBb",
                                                                 ifelse(Aa == "X" & bb == "X", "Aabb",
                                                                        ifelse(aa == "X" & BB == "X", "aaBB",
                                                                               ifelse(aa == "X" & Bb == "X", "aaBb",
                                                                                      ifelse(aa == "X" & bb == "X", "aabb", "Undefined"))))))))))

  # Define possible genotype combinations
  possible_genotypes <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")

  # Count genotypes by population
  genotype_counts <- table(Data$Pop, Data$Genotype)

  # Ensure columns are in the correct order and fill with 0 if any are missing
  genotype_counts <- as.data.frame.matrix(genotype_counts)
  for (genotype in possible_genotypes) {
    if (!genotype %in% colnames(genotype_counts)) {
      genotype_counts[[genotype]] <- 0
    }
  }

  # Add the population column
  genotype_counts$Pop <- rownames(genotype_counts)

  # Reorder the columns to match the desired order
  genotype_counts <- genotype_counts[, c("Pop", possible_genotypes)]

  # Remove row names for display purposes
  rownames(genotype_counts) <- NULL

  # Ensure all counts are numeric (convert the entire dataframe except 'Pop' column)
  genotype_counts[possible_genotypes] <- lapply(genotype_counts[possible_genotypes], as.numeric)

  # Return the dataframe
  return(genotype_counts)
}


