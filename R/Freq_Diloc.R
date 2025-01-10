#' Calcular frecuencias genotipicas y alelicas en un marcador dilocus
#'
#' Esta funcion calcula las frecuencias genotipicas y alelicas en un marcador de dos loci (dilocus),
#' dada una tabla con los genotipos y sus frecuencias observadas. Ademas, calcula las frecuencias de los haplotipos combinados
#' para los loci involucrados.
#'
#' @param tabla Un data frame que contiene dos columnas:
#' \describe{
#'   \item{Genotipo}{Vector con los genotipos observados (por ejemplo, "AABB", "AABb", "AAbb", etc.).}
#'   \item{Frecuencia}{Vector con las frecuencias observadas de cada genotipo.}
#' }
#' La tabla debe tener 9 filas correspondientes a los genotipos combinados para dos loci.
#'
#' @return Un lista con los siguientes elementos:
#' \describe{
#'   \item{Frecuencias genotipicas}{El data frame original con la columna adicional de frecuencias genotipicas.}
#'   \item{Frecuencia alelica locus A/a}{Un vector con las frecuencias alelicas para el locus A/a (p_A y q_A).}
#'   \item{Frecuencia alelica locus B/b}{Un vector con las frecuencias alelicas para el locus B/b (p_B y q_B).}
#'   \item{Frecuencia alelica Haplotipo}{Un vector con las frecuencias de los haplotipos AB, Ab, aB y ab.}
#' }
#'
#' @examples
#'
#' # Datos de entrada
#' dilocus_table <- data.frame(
#'  Pop = c("Cur1", "ANT2"),
#'  AABB = c(0, 1),
#'  AABb = c(0, 2),
#'  AAbb = c(0, 0),
#'  AaBB = c(38, 30),
#'  AaBb = c(0, 0),
#'  Aabb = c(0, 0),
#'  aaBB = c(2, 3),
#'  aaBb = c(0, 0),
#'  aabb = c(0, 0))
#'
#'  # Calcular las frecuencias dilocus para cada población
#'  resultados <- Freq_Diloc(dilocus_table)
#'
#'  # Mostrar los resultados
#'  print(resultados)
#' @export
Freq_Diloc <- function(tabla) {
  # Validación: asegurar que la tabla tenga las columnas necesarias
  RG <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb", "Pop")

  # Verificar si todas las columnas requeridas están presentes
  if (!all(RG %in% colnames(tabla))) {
    stop("La tabla debe contener las columnas de genotipos: 'AABB', 'AABb', 'AAbb', 'AaBB', 'AaBb', 'Aabb', 'aaBB', 'aaBb', 'aabb', 'Pop'")
  }

  # Inicializar lista de resultados
  resultados_lista <- list()

  # Iterar sobre cada fila (población)
  for (i in 1:nrow(tabla)) {
    fila <- tabla[i, ]

    # Calcular el total de individuos sumando las frecuencias de las columnas específicas
    total_individuos <- sum(unlist(fila[RG[-length(RG)]]))

    # Calcular las frecuencias genotípicas
    frecuencias_genotipicas <- as.numeric(fila[RG[-length(RG)]]) / total_individuos

    # Extraer las frecuencias de cada genotipo directamente desde la fila
    f_AA_BB <- as.numeric(fila["AABB"])
    f_AA_Bb <- as.numeric(fila["AABb"])
    f_AA_bb <- as.numeric(fila["AAbb"])
    f_Aa_BB <- as.numeric(fila["AaBB"])
    f_Aa_Bb <- as.numeric(fila["AaBb"])
    f_Aa_bb <- as.numeric(fila["Aabb"])
    f_aa_BB <- as.numeric(fila["aaBB"])
    f_aa_Bb <- as.numeric(fila["aaBb"])
    f_aa_bb <- as.numeric(fila["aabb"])

    # Calcular frecuencias alélicas para el locus A/a
    p_A <- (2 * (f_AA_BB + f_AA_Bb + f_AA_bb) + (f_Aa_BB + f_Aa_Bb + f_Aa_bb)) / (2 * total_individuos)
    q_A <- 1 - p_A

    # Calcular frecuencias alélicas para el locus B/b
    p_B <- (2 * (f_AA_BB + f_Aa_BB + f_aa_BB) + (f_AA_Bb + f_Aa_Bb + f_aa_Bb)) / (2 * total_individuos)
    q_B <- 1 - p_B

    # Calcular las frecuencias de haplotipos
    # Haplotipo AB
    h_AB <- (f_AA_BB + 0.5 * (f_AA_Bb + f_Aa_BB)) / total_individuos

    # Haplotipo Ab
    h_Ab <- (f_AA_bb + 0.5 * (f_AA_Bb + f_Aa_bb)) / total_individuos

    # Haplotipo aB
    h_aB <- (f_aa_BB + 0.5 * (f_Aa_BB + f_aa_Bb)) / total_individuos

    # Haplotipo ab
    h_ab <- (f_aa_bb + 0.5 * (f_Aa_bb + f_aa_Bb)) / total_individuos

    # Crear una tabla para frecuencias alélicas y de haplotipos
    freq_alelicas <- data.frame(
      "Pop" = fila["Pop"],
      "AABB" = f_AA_BB / total_individuos,
      "AABb" = f_AA_Bb / total_individuos,
      "AAbb" = f_AA_bb / total_individuos,
      "AaBB" = f_Aa_BB / total_individuos,
      "AaBb" = f_Aa_Bb / total_individuos,
      "Aabb" = f_Aa_bb / total_individuos,
      "aaBB" = f_aa_BB / total_individuos,
      "aaBb" = f_aa_Bb / total_individuos,
      "aabb" = f_aa_bb / total_individuos,
      "p_A" = p_A,
      "q_a" = q_A,
      "p_B" = p_B,
      "q_b" = q_B,
      "AB" = h_AB,
      "Ab" = h_Ab,
      "aB" = h_aB,
      "ab" = h_ab
    )

    # Añadir resultados de la fila actual a la lista de resultados
    resultados_lista[[paste("Pop", i)]] <- freq_alelicas
  }

  # Convertir lista de resultados en un data.frame
  resultados_df <- do.call(rbind, resultados_lista)
  rownames(resultados_df) <- NULL

  return(resultados_df)
}



# devtools::load_all()

# devtools::document()
