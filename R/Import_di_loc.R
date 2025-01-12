#' Import two Locus Data Mutation KDR AA/Aa/aa and BB/Bb/bb
#'
#' This function checks the consistency of a genotype matrix and counts the number of individuals with each genotype (AA, Aa, aa, BB, Bb, bb) in a dataset. It is specifically designed for working with KDR mutation genotypes and an additional locus for BB/Bb/bb genotypes, often used in resistance studies.
#'
#' @param df A data.frame called "Genotype Matrix" that represents the presence of genotypes (AA, Aa, aa, BB, Bb, bb) for a set of individuals identified by their ID and population.
#'
#' @return A data.frame with the frequency of each genotype (AA, Aa, aa, BB, Bb, bb) per population.
#'
#' @details
#' The function performs several checks to ensure data integrity:
#' - Validates that the input is a `data.frame`.
#' - Ensures the presence of the required columns: `ID`, `Pop`, `AA`, `Aa`, `aa`, `BB`, `Bb`, and `bb`.
#' - Checks for the absence of `NA` values in the required columns.
#' - Ensures that each genotype column only contains the values: `X` (indicating the individual has the genotype) and `0` (indicating absence of the genotype).
#' - Ensures that each individual (row) has exactly one `X`, meaning an individual can only present one genotype.
#' - Validates that the total count of `X` matches the total number of individuals (rows) in the data frame.
#'
#' The function counts the number of `X` for each genotype in each population and returns a frequency matrix of genotypes by population.
#'
#' @examples
#'
#' # Example 1: Provides a matrix indicating the genotype for 2 loci for 12 individuals
#' # from 1 population
#'
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
#' Import_di_loc(df)
#'
#'
#' # Example 2: Provides a matrix indicating the genotype for 106 individuals in 3 populations
#'
#' Import_di_loc(Loc2Pops3)
#'
#' @export
Import_di_loc <- function(df) {
  # Validate the data type
  if (!is.data.frame(df)) {
    stop("Error: The input must be a data.frame")
  }

  # Check column order
  correct_order <- c("ID", "Pop", "AA", "Aa", "aa", "BB", "Bb", "bb")
  if (!identical(colnames(df), correct_order)) {
    stop("Error: The columns must be in the following order: ID, Pop, AA, Aa, aa, BB, Bb, bb")
  }

  # Check for NAs in the dataframe
  if (any(is.na(df))) {
    cols_with_na <- colnames(df)[colSums(is.na(df)) > 0]
    stop(paste("Error: The data.frame contains NA values in the following columns:", paste(cols_with_na, collapse=", ")))
  }

  # Check for invalid characters in columns 3:5
  invalid_chars_1 <- df[, c("AA", "Aa", "aa")] != "X" & df[, c("AA", "Aa", "aa")] != "0"
  if (any(invalid_chars_1)) {
    cols_with_invalid_1 <- colnames(df)[3:5][apply(invalid_chars_1, 2, any)]
    stop(paste("Error: The data.frame contains invalid characters in the following columns:", paste(cols_with_invalid_1, collapse=", "), ". Only 'X' and '0' are allowed"))
  }

  # Check for invalid characters in columns 6:8
  invalid_chars_2 <- df[, c("BB", "Bb", "bb")] != "X" & df[, c("BB", "Bb", "bb")] != "0"
  if (any(invalid_chars_2)) {
    cols_with_invalid_2 <- colnames(df)[6:8][apply(invalid_chars_2, 2, any)]
    stop(paste("Error: The data.frame contains invalid characters in the following columns:", paste(cols_with_invalid_2, collapse=", "), ". Only 'X' and '0' are allowed"))
  }

  # Check for more than one "X" in each row and specify row and columns with error
  rows_with_errors <- list()
  for (i in 1:nrow(df)) {
    count_X_1 <- sum(df[i, c("AA", "Aa", "aa")] == "X") # Count the number of "X" in columns AA, Aa, aa of row i
    count_X_2 <- sum(df[i, c("BB", "Bb", "bb")] == "X") # Count the number of "X" in columns BB, Bb, bb of row i
    if (count_X_1 != 1 | count_X_2 != 1) { # If there is not exactly one "X" in either group
      rows_with_errors[[length(rows_with_errors) + 1]] <- list(row = i, count_X_1 = count_X_1, count_X_2 = count_X_2)
    }
  }

  if (length(rows_with_errors) > 0) {
    error_messages <- sapply(rows_with_errors, function(x) {
      paste("Row", x$row, "has", x$count_X_1, "'X' in columns 3:5 and", x$count_X_2, "'X' in columns 6:8")
    })
    stop(paste("Error: The following rows do not have exactly one 'X' in both groups:\n", paste(error_messages, collapse = "\n")))
  }

  # Count the total number of "X" in columns 3:5 (AA, Aa, aa)
  total_X_AA <- sum(df$AA == "X")
  total_X_Aa <- sum(df$Aa == "X")
  total_X_aa <- sum(df$aa == "X")
  total_X_group_1 <- total_X_AA + total_X_Aa + total_X_aa

  # Check if the count of 'X' in the first group (3:5) equals the number of rows
  if (total_X_group_1 != nrow(df)) {
    # Output total count of X and number of rows for debugging
    print(paste("Total X count in AA:", total_X_AA))
    print(paste("Total X count in Aa:", total_X_Aa))
    print(paste("Total X count in aa:", total_X_aa))
    print(paste("Total X count group 1 (AA, Aa, aa):", total_X_group_1))
    print(paste("Number of rows:", nrow(df)))

    stop("Error: The count of 'X' in the first group (AA, Aa, aa) must be equal to the number of rows")
  }

  # Count the total number of "X" in columns 6:8 (BB, Bb, bb)
  total_X_BB <- sum(df$BB == "X")
  total_X_Bb <- sum(df$Bb == "X")
  total_X_bb <- sum(df$bb == "X")
  total_X_group_2 <- total_X_BB + total_X_Bb + total_X_bb

  # Check if the count of 'X' in the second group (6:8) equals the number of rows
  if (total_X_group_2 != nrow(df)) {
    # Output total count of X and number of rows for debugging
    print(paste("Total X count in BB:", total_X_BB))
    print(paste("Total X count in Bb:", total_X_Bb))
    print(paste("Total X count in bb:", total_X_bb))
    print(paste("Total X count group 2 (BB, Bb, bb):", total_X_group_2))
    print(paste("Number of rows:", nrow(df)))

    stop("Error: The count of 'X' in the second group (BB, Bb, bb) must be equal to the number of rows")
  }

  # Create an empty data frame to store frequencies by population
  freq_matrix <- data.frame(
    Pop = character(length(unique(df$Pop))),
    AA = numeric(length(unique(df$Pop))),
    Aa = numeric(length(unique(df$Pop))),
    aa = numeric(length(unique(df$Pop))),
    BB = numeric(length(unique(df$Pop))),
    Bb = numeric(length(unique(df$Pop))),
    bb = numeric(length(unique(df$Pop))),
    stringsAsFactors = FALSE
  )

  # Get unique populations
  populations <- unique(df$Pop)

  for (i in seq_along(populations)) {
    pop <- populations[i]
    df_pop <- df[df$Pop == pop, ]

    # Count the number of individuals with each genotype for the population
    AA <- sum(df_pop$AA == "X")
    Aa <- sum(df_pop$Aa == "X")
    aa <- sum(df_pop$aa == "X")
    BB <- sum(df_pop$aa == "X")
    Bb <- sum(df_pop$Bb == "X")
    bb <- sum(df_pop$bb == "X")

    # Store the results in the data frame
    freq_matrix[i, ] <- c(pop, AA, Aa, aa, BB, Bb, bb)
  }


  # If all checks pass, proceed with the second analysis using GenDiLoc
  result <- GenDiLoc(df)

  return(list(MolocusGenotipes = data.frame(freq_matrix), DilocusGenotypesMatrix = data.frame(result)))
}
