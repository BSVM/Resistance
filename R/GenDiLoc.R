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

  # Check values are only '0' or 'X' in genotype columns
  genotype_cols <- setdiff(names(Data), c("ID", "Pop"))
  valid_values <- c("0", "X")
  for (col in genotype_cols) {
    if (!all(Data[[col]] %in% valid_values)) {
      stop(paste("Column", col, "contains invalid values."))
    }
  }

  # Check each row has exactly one 'X' per locus
  Data$A_Locus <- apply(Data[, c("AA", "Aa", "aa")], 1, function(x) sum(x == "X"))
  Data$B_Locus <- apply(Data[, c("BB", "Bb", "bb")], 1, function(x) sum(x == "X"))

  if (any(Data$A_Locus != 1) || any(Data$B_Locus != 1)) {
    stop("Each row must have exactly one 'X' per locus (A and B).")
  }

  Data <- Data[, !names(Data) %in% c("A_Locus", "B_Locus")]  # Remove helper columns

  # Create Genotype column using vectorized approach
  A <- names(Data)[3:5][max.col(Data[, 3:5] == "X", "first")]
  B <- names(Data)[6:8][max.col(Data[, 6:8] == "X", "first")]

  Data$Genotype <- paste0(
    ifelse(A == "AA", "AA", ifelse(A == "Aa", "Aa", "aa")),
    ifelse(B == "BB", "BB", ifelse(B == "Bb", "Bb", "bb"))
  )

  # Ensure correct genotype names (e.g., AABb instead of AABb)
  Data$Genotype <- sub("AAbb", "AAbb", Data$Genotype)  # Fix any case inconsistencies if needed

  # Count genotypes by population
  possible_genotypes <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")
  genotype_counts <- as.data.frame.matrix(table(Data$Pop, Data$Genotype))

  # Add missing genotype columns with 0s
  missing <- setdiff(possible_genotypes, names(genotype_counts))
  genotype_counts[missing] <- 0
  genotype_counts <- genotype_counts[, possible_genotypes]  # Reorder columns

  genotype_counts$Pop <- rownames(genotype_counts)
  rownames(genotype_counts) <- NULL
  genotype_counts <- genotype_counts[, c("Pop", possible_genotypes)]

  return(genotype_counts)
}

