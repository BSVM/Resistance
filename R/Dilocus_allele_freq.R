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
#' #Crear una tabla de ejemplo con genotipos y frecuencias
#' dilocus_tabla <- data.frame(
#'   Genotipo = c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb"),
#'   Frecuencia = c(0, 0, 0, 38, 0, 0, 2, 0, 0))
#'
#
#' calcular_frecuencias_dilocus(dilocus_tabla)
#'
#' @export
calcular_frecuencias_dilocus <- function(tabla) {
  # Validación: asegurar que la tabla tenga las columnas necesarias
  if (!all(c("Genotipo", "Frecuencia") %in% colnames(tabla))) {
    stop("La tabla debe contener las columnas 'Genotipo' y 'Frecuencia'.")
  }

  # Calcular el total de individuos
  total_individuos <- sum(tabla$Frecuencia)

  # Calcular las frecuencias genotípicas
  tabla$Frecuencia_Genotipica <- tabla$Frecuencia / total_individuos

  # Extraer las frecuencias de cada genotipo
  f_AA_BB <- tabla$Frecuencia[tabla$Genotipo == "AABB"]
  f_AA_Bb <- tabla$Frecuencia[tabla$Genotipo == "AABb"]
  f_AA_bb <- tabla$Frecuencia[tabla$Genotipo == "AAbb"]
  f_Aa_BB <- tabla$Frecuencia[tabla$Genotipo == "AaBB"]
  f_Aa_Bb <- tabla$Frecuencia[tabla$Genotipo == "AaBb"]
  f_Aa_bb <- tabla$Frecuencia[tabla$Genotipo == "Aabb"]
  f_aa_BB <- tabla$Frecuencia[tabla$Genotipo == "aaBB"]
  f_aa_Bb <- tabla$Frecuencia[tabla$Genotipo == "aaBb"]
  f_aa_bb <- tabla$Frecuencia[tabla$Genotipo == "aabb"]

  # Calcular frecuencias alelicas para el locus A/a
  p_A <- (2 * (f_AA_BB + f_AA_Bb + f_AA_bb) + (f_Aa_BB + f_Aa_Bb + f_Aa_bb)) / (2 * total_individuos)
  q_A <- 1 - p_A

  # Calcular frecuencias alelicas para el locus B/b
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

  # Resultados detallados
  resultados <- list(
    "Frecuencias genotipicas" = tabla,
    "Frecuencia alelica locus A/a" = list("p_A (A)" = p_A, "q_A (a)" = q_A),
    "Frecuencia alelica locus B/b" = list("p_B (B)" = p_B, "q_B (b)" = q_B),
    "Frecuencia alelica Haplotype" = list("V/F" = h_ab, "V/C" = h_aB, "I/C" = h_AB, "I/F" = h_Ab)
  )

  return(resultados)
}


