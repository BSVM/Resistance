#' Calcular frecuencias genotipicas y alelicas sin locus
#'
#' Esta funcion calcula las frecuencias genotipicas y alelicas a partir de un vector que contiene los conteos de individuos para cada genotipo (AA, Aa, aa).
#'
#' @param genotype_counts Un vector numerico de longitud 3 que representa los conteos de individuos con los genotipos: AA, Aa y aa, en ese orden.
#'
#' @return Una lista con tres elementos:
#' \describe{
#'   \item{genotypic_frequencies}{Un data.frame con las frecuencias genotipicas de AA, Aa, y aa.}
#'   \item{allelic_frequencies}{Un data.frame con las frecuencias alelicas de A y a.}
#'   \item{total_individuals}{El numero total de individuos considerados en los calculos.}
#' }
#'
#' @examples
#'
#' genotype_counts <- c(AA=50, Aa=30, aa=20)  # 50 individuos con genotipo AA, 30 con Aa, y 20 con aa
#'
#' result <- Single_locus_frequencies(genotype_counts)
#'
#' print(result$genotypic_frequencies) # Frecuencias genotipicas
#' print(result$allelic_frequencies)  # Frecuencias alelicas
#' print(result$total_individuals)    # Total de individuos
#'
#' @export
Single_locus_frequencies <- function(genotype_counts) {
  # Extraer conteos de cada genotipo
  f_AA <- genotype_counts[1]
  f_Aa <- genotype_counts[2]
  f_aa <- genotype_counts[3]

  # Calcular el numero total de individuos
  total_individuals <- sum(genotype_counts)

  # Calcular las frecuencias genotipicas
  freq_AA <- f_AA / total_individuals
  freq_Aa <- f_Aa / total_individuals
  freq_aa <- f_aa / total_individuals

  # Calcular las frecuencias alelicas
  p <- (f_AA + 0.5 * f_Aa) / total_individuals  # Frecuencia del alelo A
  q <- 1 - p  # Frecuencia del alelo a

  # Devolver las frecuencias y el total de individuos
  list(
    genotypic_frequencies = data.frame(
      Genotype = c("AA", "Aa", "aa"),
      Frequency = c(freq_AA, freq_Aa, freq_aa)
    ),
    allelic_frequencies = data.frame(
      Allele = c("A", "a"),
      Frequency = c(p, q)
    ),
    total_individuals = total_individuals
  )
}

# devtools::load_all()

# devtools::document()
