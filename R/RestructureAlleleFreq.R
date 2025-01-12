#' Restructure Allele Frequency Data
#'
#' This function restructures allele frequency data into a long format. It accepts data in two formats, specified by the `Loc` parameter.
#' In `Loc = 1`, the input data should have columns for the frequencies of two alleles (Freq_A and Freq_a). In `Loc = 2`, the data
#' should include columns for the frequencies of four genotypes (AB, Ab, aB, ab).
#'
#' @param data A data frame containing population data and allele frequencies. The structure of the data frame depends on the value of the `Loc` parameter.
#' @param Loc An integer indicating the format of the data:
#' 1: The data frame should have columns "Pop", "Freq_A", and "Freq_a" for allele frequencies.
#' 2: The data frame should have columns "Pop", "AB", "Ab", "aB", and "ab" for genotype frequencies. Default is 1.
#'
#' @return A data frame in long format with three columns: `Population`, `Allele`, and `Frequency`.
#' The function returns a long-format data frame where each row represents an allele and its corresponding frequency in a population.
#'
#' @export
#'
#' @examples
#' # Example 1 for Loc = 1
#' df1 <- data.frame(
#'   Pop = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6"),
#'   Freq_A = c(0.6, 0.7, 0.5, 0.8, 0.4, 0.3),
#'   Freq_a = c(0.4, 0.3, 0.5, 0.2, 0.6, 0.7))
#'
#' allele_data_1 <- RestructureAlleleFreq(df1, Loc = 1)
#'
#' # Example 2 for Loc = 2
#' df2 <- data.frame(
#'   Pop = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6"),
#'   AB = c(0.2, 0.3, 0.1, 0.25, 0.4, 0.35),
#'   Ab = c(0.3, 0.4, 0.2, 0.15, 0.25, 0.3),
#'   aB = c(0.4, 0.2, 0.3, 0.45, 0.2, 0.25),
#'   ab = c(0.1, 0.1, 0.4, 0.15, 0.15, 0.1))
#'
#' allele_data_2 <- RestructureAlleleFreq(df2, Loc = 2)
RestructureAlleleFreq <- function(data, Loc = 1) {
  # Data input checks
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }

  if (Loc == 1) {
    required_columns <- c("Pop", "Freq_A", "Freq_a")
    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      stop(paste("The data frame must contain the following columns:", paste(missing_columns, collapse = ", ")))
    }

    # Check that the frequency columns are numeric
    if (!all(sapply(data[, c("Freq_A", "Freq_a")], is.numeric))) {
      stop("The columns 'Freq_A' and 'Freq_a' must contain numeric values.")
    }

    # Check that the frequencies are between 0 and 1
    if (any(data$Freq_A < 0 | data$Freq_A > 1 | data$Freq_a < 0 | data$Freq_a > 1)) {
      stop("Allele frequencies must be between 0 and 1.")
    }

    # Check that allele frequencies sum to 1 for each population
    if (any(abs(data$Freq_A + data$Freq_a - 1) > 1e-6)) {
      stop("The frequencies of alleles A and a must sum to 1 for each population.")
    }

    # Check that the Pop column has no missing values
    if (any(is.na(data$Pop))) {
      stop("The 'Pop' column must not contain missing values.")
    }

    # Restructure the data into long format
    allele_data <- data.frame(
      Population = rep(data$Pop, each = 2),
      Allele = rep(c("A", "a"), times = nrow(data)),
      Frequency = c(t(data[, c("Freq_A", "Freq_a")]))  # Accessing the frequencies correctly
    )

  } else if (Loc == 2) {
    required_columns <- c("Pop", "AB", "Ab", "aB", "ab")
    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      stop(paste("The data frame must contain the following columns:", paste(missing_columns, collapse = ", ")))
    }

    # Check that the frequency columns are numeric
    if (!all(sapply(data[, c("AB", "Ab", "aB", "ab")], is.numeric))) {
      stop("The columns 'AB', 'Ab', 'aB', and 'ab' must contain numeric values.")
    }

    # Check that the frequencies are between 0 and 1
    if (any(data$AB < 0 | data$AB > 1 | data$Ab < 0 | data$Ab > 1 | data$aB < 0 | data$aB > 1 | data$ab < 0 | data$ab > 1)) {
      stop("Allele frequencies must be between 0 and 1.")
    }

    # Check that the Pop column has no missing values
    if (any(is.na(data$Pop))) {
      stop("The 'Pop' column must not contain missing values.")
    }

    # Restructure the data into long format
    allele_data <- data.frame(
      Population = rep(data$Pop, each = 4),
      Allele = rep(c("AB", "Ab", "aB", "ab"), times = nrow(data)),
      Frequency = c(t(data[, c("AB", "Ab", "aB", "ab")]))  # Accessing the frequencies correctly
    )

  } else {
    stop("Invalid value for 'Loc'. It must be 1 or 2.")
  }

  return(allele_data)
}

