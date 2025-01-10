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
#' # Example 1: Provides a matrix indicating the genotype for 12 individuals
#'
#' df <- data.frame(
#'   ID = 1:12,
#'   aa = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
#'   Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
#'   AA = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Import_single_loc(df)
#'
#' @export
Import_single_loc <- function(df) {
  # Validate data type
  if (!is.data.frame(df)) {
    stop("Error: The input must be a data.frame")
  }

  # Check required columns
  required_cols <- c("aa", "Aa", "AA")
  if (!all(required_cols %in% colnames(df))) {
    stop("Error: The data.frame must contain the columns 'aa', 'Aa', and 'AA'")
  }

  # Check for NAs in the dataframe
  if (any(is.na(df[required_cols]))) {
    stop("Error: The data.frame contains NA values")
  }

  # Check for invalid characters
  invalid_chars <- df[required_cols] != "X" & df[required_cols] != "0"
  if (any(invalid_chars)) {
    stop("Error: The data.frame contains invalid characters. Only 'X' and '0' are allowed")
  }

  # Check if there is more than one X in each row
  for (i in 1:nrow(df)) {
    if (sum(df[i, required_cols] == "X") > 1) {
      stop(paste("Error: There is more than one 'X' in row", i))
    }
  }

  # Count the total number of "X"
  total_X <- sum(df[required_cols] == "X")

  # Check if the count of "X" is equal to the number of rows
  if (total_X != nrow(df)) {
    stop("Error: The count of 'X' must be equal to the number of rows")
  }

  # Count the number of individuals with each genotype
  aa <- sum(df$aa == "X")
  Aa <- sum(df$Aa == "X")
  AA <- sum(df$AA == "X")

  # Create the frequency table
  freq_table <- as.numeric(c(AA, Aa, aa))

  return(freq_table)
}

