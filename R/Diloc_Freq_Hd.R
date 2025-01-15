#' Calculate Allele Frequencies for a Biallelic Marker Across Multiple Populations
#'
#' This function computes the haplotype frequencies for a biallelic genetic marker (AB, Ab, aB, ab) in multiple populations
#' based on a genotype distribution. The method accounts for genotype counts and utilizes an iterative approach to
#' estimate the haplotype frequencies. The results can be returned as the final iteration or with the full iteration history
#' for each population. This function is designed to work with a data frame that contains genotype counts for each population.
#'
#' @param data A data frame where each row corresponds to a different population. The first column should contain the population names
#' (`Pop`), and the remaining columns should represent the genotype counts for the following genotypes: `AABB`, `AABb`, `AAbb`, `AaBB`,
#' `AaBb`, `Aabb`, `aaBB`, `aaBb`, `aabb`.
#'
#' @param max_iterations The maximum number of iterations to perform in the estimation process. Default is 15.
#'
#' @param tol The tolerance level for convergence. If the change in the correction factors (`CorAB_ab` and `CorAb_aB`)
#' between iterations is smaller than this value, the iteration will stop. Default is 1e-6.
#'
#' @param dat Specifies which iteration results to return:
#'   - `dat = 1` returns only the final results of the last iteration for each population.
#'   - `dat = 2` returns the full iteration history for each population.
#'   Default is `dat = 1`.
#'
#' @return A data frame containing the haplotype frequencies (`AB`, `Ab`, `aB`, `ab`) for each population.
#' The output also includes the iteration number, population name, total sample size (`N`), and the correction factors
#' (`CorAB_ab`, `CorAb_aB`) for each iteration.
#'
#' @details The function uses a correction approach to estimate haplotype frequencies.
#' The iterative process aims to converge by adjusting the correction factors (`CorAB_ab` and `CorAb_aB`) until their change
#' between iterations is smaller than the specified tolerance (`tol`). The iteration history can be useful for assessing convergence
#' behavior or troubleshooting the process.
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
#' Diloc_Freq_Hd(DicGtps, max_iterations = 10, tol = 1e-6, dat = 1)
#'
#' @export
Diloc_Freq_Hd <- function(data, max_iterations = 10, tol = 1e-6, dat = 1) {
  # Check if the "Pop" column exists
  if (!"Pop" %in% colnames(data)) {
    stop("The 'Pop' column does not exist in the data frame.")
  }

  # Check if there are any NA values in the data frame
  if (any(is.na(data))) {
    stop("The data frame contains NA or missing values.")
  }

  # Check if the necessary columns are present
  required_columns <- c("AABB", "AABb", "AAbb", "AaBB", "AaBb", "Aabb", "aaBB", "aaBb", "aabb")
  if (!all(required_columns %in% colnames(data))) {
    stop("The data frame does not contain all the necessary columns.")
  }

  # Check if the genotype columns are numeric and non-negative
  if (any(sapply(data[required_columns], function(x) !is.numeric(x) || any(x < 0)))) {
    stop("The genotype columns must be numeric and contain no negative values.")
  }

  # Initialize a list to store the results of all iterations
  all_results <- list()

  # Identify the genotype columns
  N_c <- required_columns

  # Define an internal function to calculate haplotype frequencies
  calc_freqs <- function(data, i, CorAB_ab, CorAb_aB) {
    N <- sum(data[i, N_c])
    f_AB <- (2 * data$AABB[i] + data$AABb[i] + data$AaBB[i] + CorAB_ab * data$AaBb[i]) / (2 * N)
    f_Ab <- (2 * data$AAbb[i] + data$AABb[i] + data$Aabb[i] + CorAb_aB * data$AaBb[i]) / (2 * N)
    f_aB <- (2 * data$aaBB[i] + data$AaBB[i] + data$aaBb[i] + CorAb_aB * data$AaBb[i]) / (2 * N)
    f_ab <- (2 * data$aabb[i] + data$Aabb[i] + data$aaBb[i] + CorAB_ab * data$AaBb[i]) / (2 * N)
    return(c(f_AB, f_Ab, f_aB, f_ab))
  }

  # Iterate over each row (population)
  for (i in 1:nrow(data)) {
    # Initialize correction values for the first iteration
    CorAB_ab <- 0.5
    CorAb_aB <- 0.5

    iter_results <- data.frame(Iter = integer(), Pop = character(), N = integer(), AB = numeric(), Ab = numeric(), aB = numeric(), ab = numeric(), CorAB_ab = numeric(), CorAb_aB = numeric())

    # Save the initial population values
    initial_freqs <- calc_freqs(data, i, CorAB_ab, CorAb_aB)

    for (iter in 1:max_iterations) {
      N <- sum(data[i, N_c])

      # Calculate the haplotype frequencies
      freqs <- calc_freqs(data, i, CorAB_ab, CorAb_aB)

      # Calculate the ambiguity error for each population independently
      k1 <- 2 * freqs[1] * freqs[4] # AB/ab
      k2 <- 2 * freqs[2] * freqs[3] # Ab/aB

      # Check if k1 + k2 is zero, which may be problematic
      if ((k1 + k2) != 0) {
        new_CorAB_ab <- k1 / (k1 + k2) # AB/ab
        new_CorAb_aB <- k2 / (k1 + k2) # Ab/aB
      } else {
        new_CorAB_ab <- NA
        new_CorAb_aB <- NA
        # Display initial frequency values for the population with the error
        message(paste("Population", data$Pop[i], "has problematic frequency values in iteration", iter))
        message(paste("Initial frequencies: AB =", initial_freqs[1], "Ab =", initial_freqs[2], "aB =", initial_freqs[3], "ab =", initial_freqs[4]))

        # Store the initial population values in the final results
        iter_results <- rbind(iter_results, data.frame(Iter = "Initial", Pop = data$Pop[i], N = N, AB = initial_freqs[1], Ab = initial_freqs[2], aB = initial_freqs[3], ab = initial_freqs[4], CorAB_ab = NA, CorAb_aB = NA))
        break  # Exit the loop and continue with the next population
      }

      # Store the results of the current iteration
      iter_results <- rbind(iter_results, data.frame(Iter = iter, Pop = data$Pop[i], N = N, AB = freqs[1], Ab = freqs[2], aB = freqs[3], ab = freqs[4], CorAB_ab = new_CorAB_ab, CorAb_aB = new_CorAb_aB))

      # Check for convergence
      if (!is.na(new_CorAB_ab) && !is.na(new_CorAb_aB) && abs(new_CorAB_ab - CorAB_ab) < tol && abs(new_CorAb_aB - CorAb_aB) < tol) {
        break
      }

      # Update correction values for the next iteration
      CorAB_ab <- new_CorAB_ab
      CorAb_aB <- new_CorAb_aB
    }

    # Store the results of all iterations or only the last one, depending on the `dat` parameter
    if (dat == 2) {
      all_results[[data$Pop[i]]] <- iter_results
    } else {
      all_results[[data$Pop[i]]] <- iter_results[nrow(iter_results), ]
    }
  }

  # Combine the results of all populations into a single data frame
  final_results <- do.call(rbind, all_results)

  return(final_results)
}
