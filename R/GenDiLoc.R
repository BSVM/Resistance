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
  # Create a new column 'Genotype' based on the conditions
  Data$Genotype <- with(Data, ifelse(AA == "X" & BB == "X", "AABB",
                                     ifelse(AA == "X" & Bb == "X", "AABb",
                                            ifelse(AA == "X" & bb == "X", "AAbb",
                                                   ifelse(Aa == "X" & BB == "X", "AaBB",
                                                          ifelse(Aa == "X" & Bb == "X", "AaBb",
                                                                 ifelse(Aa == "X" & bb == "X", "Aabb",
                                                                        ifelse(aa == "X" & BB == "X", "aaBB",
                                                                               ifelse(aa == "X" & Bb == "X", "aaBb",
                                                                                      ifelse(aa == "X" & bb == "X", "aabb", "Indefinido"))))))))))

  # Define the possible genotype combinations
  genotipos_posibles <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")

  # Count genotypes by population
  conteo_genotipos <- table(Data$Pop, Data$Genotype)

  # Ensure columns are in the correct order and fill with 0 if any are missing
  conteo_genotipos <- as.data.frame.matrix(conteo_genotipos)
  for (genotipo in genotipos_posibles) {
    if (!genotipo %in% colnames(conteo_genotipos)) {
      conteo_genotipos[[genotipo]] <- 0
    }
  }

  # Add the population column
  conteo_genotipos$Pop <- rownames(conteo_genotipos)

  # Reorder the columns to match the desired order
  conteo_genotipos <- conteo_genotipos[, c("Pop", genotipos_posibles)]

  # Remove row names for display purposes
  rownames(conteo_genotipos) <- NULL

  # Convert the data frame to a matrix
  conteo_genotipos_matrix <- as.matrix(conteo_genotipos)

  # Return the matrix
  return(conteo_genotipos_matrix)
}
