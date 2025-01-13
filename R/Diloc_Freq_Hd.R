#' Calculate Allele Frequencies for a Biallelic Marker
#'
#' This function calculates haplotype frequencies for a biallelic genetic marker based on a genotype distribution
#' for multiple populations. The input is a data frame with the genotype counts for each population.
#'
#' @param data A data frame where each row represents a population, the first column is the population name (`Pop`),
#' and the remaining columns are genotype counts (`AABB`, `AABb`, `AAbb`, `AaBB`, `AaBb`, `Aabb`, `aaBB`, `aaBb`, `aabb`).
#'
#' @return A data frame with haplotype frequencies (`AB`, `Ab`, `aB`, `ab`) for each population.
#'
#' @examples
#' # Create a data frame with genotype counts for three populations
#' DicGtps <- data.frame(
#'   Pop = c("Pop 1", "Pop 2", "Pop 3"),
#'   AABB = c(5, 3, 4),
#'   AABb = c(3, 2, 5),
#'   AAbb = c(4, 2, 5),
#'   AaBB = c(4, 6, 3),
#'   AaBb = c(10, 8, 12),
#'   Aabb = c(6, 5, 7),
#'   aaBB = c(3, 4, 2),
#'   aaBb = c(2, 1, 3),
#'   aabb = c(3, 2, 4)
#' )
#'
#' # Calculate haplotype frequencies for the data frame
#' Diloc_Freq_Hd(DicGtps)
#'
#' @export
Diloc_Freq_Hd <- function(data) {
  # Ensure the "Pop" column exists
  if (!"Pop" %in% colnames(data)) {
    stop("The 'Pop' column is missing from the data frame.")
  }

  # Check for missing or NA values
  if (any(is.na(data))) {
    stop("The data frame contains missing (NA) values.")
  }

  # Check for required genotype columns
  required_columns <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")
  if (!all(required_columns %in% colnames(data))) {
    stop("The data frame must contain all genotype columns: ", paste(required_columns, collapse = ", "))
  }

  # Ensure genotype columns are numeric and non-negative
  if (any(sapply(data[required_columns], function(x) !is.numeric(x) || any(x < 0)))) {
    stop("Genotype columns must be numeric and non-negative.")
  }

  # Initialize a list to store results
  results <- list()

  # Loop through each population to calculate haplotype frequencies
  for (i in seq_len(nrow(data))) {
    # Total number of individuals in the population
    total_individuals <- sum(data[i, required_columns])

    # Check for zero total individuals
    if (total_individuals == 0) {
      stop(paste("Total individuals in population", data$Pop[i], "is zero."))
    }

    # Calculate haplotype frequencies
    f_AB <- (2 * data$AABB[i] + data$AABb[i] + data$AaBB[i] + data$AaBb[i]) / (2 * total_individuals)
    f_Ab <- (2 * data$AAbb[i] + data$AABb[i] + data$Aabb[i] + data$AaBb[i]) / (2 * total_individuals)
    f_aB <- (2 * data$aaBB[i] + data$AaBB[i] + data$aaBb[i] + data$AaBb[i]) / (2 * total_individuals)
    f_ab <- (2 * data$aabb[i] + data$Aabb[i] + data$aaBb[i] + data$AaBb[i]) / (2 * total_individuals)

    # Store results
    results[[data$Pop[i]]] <- c(f_AB, f_Ab, f_aB, f_ab)
  }

  # Convert results to a data frame
  return(data.frame(
    Pop = data$Pop,
    AB = sapply(results, `[[`, 1),
    Ab = sapply(results, `[[`, 2),
    aB = sapply(results, `[[`, 3),
    ab = sapply(results, `[[`, 4)
  ))
}
