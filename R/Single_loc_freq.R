#' Calculate Genotypic and Allelic Frequencies
#'
#' This function calculates the genotypic and allelic frequencies from a matrix containing genotype counts (AA, Aa, aa) for each population.
#'
#' @param freq_matrix A data frame with the following columns:
#' \describe{
#'   \item{Pop}{A unique identifier for each population (e.g., "Pop1", "Pop2").}
#'   \item{aa}{The count of individuals with the genotype aa.}
#'   \item{Aa}{The count of individuals with the genotype Aa.}
#'   \item{AA}{The count of individuals with the genotype AA.}
#' }
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Population}{The identifier for each population.}
#'   \item{n}{The total number of individuals considered in the calculations.}
#'   \item{aa}{The genotypic frequency of individuals with the aa genotype.}
#'   \item{Aa}{The genotypic frequency of individuals with the Aa genotype.}
#'   \item{AA}{The genotypic frequency of individuals with the AA genotype.}
#'   \item{Freq_A}{The allelic frequency of allele A.}
#'   \item{Freq_a}{The allelic frequency of allele a.}
#' }
#'
#' @examples
#' # Example: Genotypic frequency matrix for a single population
#' freq_matrix <- data.frame(
#'   Pop = c("Pop1"),
#'   aa = c(2),
#'   Aa = c(9),
#'   AA = c(1))
#'
#' # Calculate frequencies
#' Single_loc_freq(freq_matrix)
#'
#' # Example: Genotypic frequency matrix for multiple populations
#' df2 <- Import_single_loc(Loc1Pops3)
#' Single_loc_freq(df2)
#'
#' @export
Single_loc_freq <- function(freq_matrix) {
  # Check if the input is a data frame
  if (!is.data.frame(freq_matrix)) {
    stop("The argument 'freq_matrix' must be a data frame.")
  }

  # Required columns
  required_columns <- c("Pop", "aa", "Aa", "AA")
  if (!all(required_columns %in% colnames(freq_matrix))) {
    stop("The frequency matrix must contain the columns: 'Pop', 'aa', 'Aa', and 'AA'.")
  }

  # Check data types
  if (!all(sapply(freq_matrix[, c("aa", "Aa", "AA")], function(x) is.numeric(x) && all(x == floor(x))))) {
    stop("The values in the columns 'aa', 'Aa', and 'AA' must be integers.")
  }

  # Additional validations
  if (any(freq_matrix[, c("AA", "Aa", "aa")] < 0)) {
    stop("The values in the columns 'AA', 'Aa', and 'aa' cannot be negative.")
  }
  if (any(is.na(freq_matrix[, c("AA", "Aa", "aa")]))) {
    stop("The values in the columns 'AA', 'Aa', and 'aa' cannot be NA.")
  }
  if (any(rowSums(freq_matrix[, c("AA", "Aa", "aa")]) == 0)) {
    stop("The total number of individuals per population cannot be zero.")
  }
  if (any(is.na(freq_matrix$Pop) | freq_matrix$Pop == "")) {
    stop("The 'Pop' column cannot contain empty or NA values.")
  }
  if (length(unique(freq_matrix$Pop)) != nrow(freq_matrix)) {
    stop("The population identifiers in the 'Pop' column must be unique.")
  }

  # Create the result data frame
  result_df <- data.frame(
    Pop = character(0),
    n = integer(0),
    aa = numeric(0),
    Aa = numeric(0),
    AA = numeric(0),
    Freq_A = numeric(0),
    Freq_a = numeric(0),
    stringsAsFactors = FALSE
  )

  # Calculate frequencies by population
  for (i in 1:nrow(freq_matrix)) {
    genotype_counts <- as.numeric(freq_matrix[i, c("AA", "Aa", "aa")])
    total_individuals <- sum(genotype_counts)
    freq_AA <- genotype_counts[1] / total_individuals
    freq_Aa <- genotype_counts[2] / total_individuals
    freq_aa <- genotype_counts[3] / total_individuals
    p <- (genotype_counts[1] + 0.5 * genotype_counts[2]) / total_individuals
    q <- 1 - p

    result_df <- rbind(result_df, data.frame(
      Pop = freq_matrix$Pop[i],
      n = total_individuals,
      AA = freq_AA,
      Aa = freq_Aa,
      aa = freq_aa,
      Freq_A = p,
      Freq_a = q
    ))
  }

  return(result_df)
}
