#' Plot Allele Frequencies by Population
#'
#' This function creates a bar plot of allele frequencies by population. It can handle two types of input data structures,
#' specified by the `Loc` argument.
#'
#' @param data A data frame containing the population, allele, and frequency data.
#' @param Loc An integer (1 or 2) indicating the type of input data structure.
#' \itemize{
#'   \item{1: Data structure with alleles "A" and "a".}
#'   \item{2: Data structure with haplotypes "Ab", "AB", "aB", and "ab".}
#' }
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#'
#' @return A ggplot object representing the bar plot of allele frequencies by population.
#' @export
#'
#' @examples
#'
#' # Example 1 for Loc = 1
#' data1 <- data.frame(
#'   Pop = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6"),
#'   Freq_A = c(0.6, 0.7, 0.5, 0.8, 0.4, 0.3),
#'   Freq_a = c(0.4, 0.3, 0.5, 0.2, 0.6, 0.7)
#' )
#'
#' dataf1 <- RestructureAlleleFreq(data1, Loc = 1)
#' PlotAlleleFrequency(dataf1, Loc = 1)
#'
#' # Example 2 for Loc = 2
#' data2 <- data.frame(
#'   Pop = c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6"),
#'   AB = c(0.2, 0.3, 0.1, 0.25, 0.4, 0.35),
#'   Ab = c(0.3, 0.4, 0.2, 0.15, 0.25, 0.3),
#'   aB = c(0.4, 0.2, 0.3, 0.45, 0.2, 0.25),
#'   ab = c(0.1, 0.1, 0.4, 0.15, 0.15, 0.1)
#' )
#'
#' dataf2 <- RestructureAlleleFreq(data2, Loc = 2)
#'
#' PlotAlleleFrequency(dataf2, Loc = 2)
PlotAlleleFrequency <- function(data, Loc = 1) {

  # Data input checks common for both types
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }

  if (Loc == 1) {
    # Define required columns for Loc = 1
    required_columns <- c("Population", "Allele", "Frequency")

    # Check for missing columns
    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      stop(paste("The data frame must contain the following columns:", paste(missing_columns, collapse = ", ")))
    }

    # Ensure allele levels are in the specified order for Loc = 1
    allele_order <- c("A", "a")
    data$Allele <- factor(data$Allele, levels = allele_order)

    # Create the plot for Loc = 1
    plot <- ggplot(data, aes(x = Population, y = Frequency, fill = Allele)) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.5) +
      scale_fill_manual(values = c("A" = "#FF6347", "a" = "#4682B4")) +
      labs(x = "Population", y = "Allele Frequency", fill = "Allele") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")
      )

  } else if (Loc == 2) {
    # Define required columns for Loc = 2
    required_columns <- c("Population", "Allele", "Frequency")

    # Check for missing columns
    missing_columns <- setdiff(required_columns, colnames(data))
    if (length(missing_columns) > 0) {
      stop(paste("The data frame must contain the following columns:", paste(missing_columns, collapse = ", ")))
    }

    # Ensure allele levels are in the specified order for Loc = 2
    allele_order <- c("Ab", "AB", "aB", "ab")
    data$Allele <- factor(data$Allele, levels = allele_order)

    # Create the plot for Loc = 2
    plot <- ggplot(data, aes(x = Population, y = Frequency, fill = Allele)) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.5) +
      scale_fill_manual(values = c("ab" = "#1E90FF", "aB" = "yellow", "AB" = "red", "Ab" = "darkblue")) +
      labs(x = "Population", y = "Frequency Haplotypes", fill = "Allele") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 12, face = "plain", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 12, face = "plain", angle = 0, hjust = 1),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", angle = 90, hjust = 0.5)
      )
  } else {
    stop("Invalid value for 'Loc'. It must be 1 or 2.")
  }

  # Print the plot
  print(plot)
}
